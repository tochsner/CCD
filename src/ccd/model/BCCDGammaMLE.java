package ccd.model;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;


public class BCCDGammaMLE extends ParameterEstimator<BCCD> {
    Function<CladePartitionObservation, Double>[] getBetaObservations;
    Function<CladePartitionObservation, Double> feature;
    BCCDFeatureSelector featureSelector;
    int[][] betaGroups;

    public BCCDGammaMLE(Function<CladePartitionObservation, Double> feature) {
        this.feature = feature;
    }

    public BCCDGammaMLE(BCCDFeatureSelector featureSelector) {
        this.featureSelector = featureSelector;
    }

    @Override
    public BCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new BCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(BCCD bccd) {
        List<BCCDCladePartition> partitions = bccd.getAllPartitions();

        if (this.feature != null) {
            this.getBetaObservations = new Function[partitions.size()];
            for (int i = 0; i < partitions.size(); i++) {
                this.getBetaObservations[i] = feature;
            }
        } else {
            this.getBetaObservations = this.featureSelector.getFeatures(partitions);
        }

        this.betaGroups = new int[partitions.size()][];
        for (int i = 0; i < partitions.size(); i++) {
            this.betaGroups[i] = new int[]{i};
        }

        double approximateShape = this.estimateMedianApproximateShape(partitions);

        for (int[] betaGroup : this.betaGroups) {
            BrentOptimizer solver = new BrentOptimizer(1e-3, 1e-3);
            BetaOptimizer optimizer = new BetaOptimizer(
                    bccd,
                    betaGroup,
                    this.getBetaObservations,
                    approximateShape
            );

            double optimalBeta;
            try {
                optimalBeta = solver.optimize(
                        GoalType.MAXIMIZE,
                        new InitialGuess(new double[]{0.0}),
                        new SearchInterval(-3, 3),
                        new UnivariateObjectiveFunction(optimizer),
                        new MaxEval(100)
                ).getPoint();
            } catch (NoBracketingException e) {
                throw new ArithmeticException("No solution found");
            }

            double[] optimalShapes = optimizer.getOptimalShapes(optimalBeta);
            double[] optimalScales = optimizer.getOptimalScales(optimalBeta, optimalShapes);

            for (int i = 0; i < betaGroup.length; i++) {
                setDistribution(
                        bccd, betaGroup[i], optimalBeta, optimalShapes[i], optimalScales[i]
                );
            }
        }
    }

    private double estimateMedianApproximateShape(List<BCCDCladePartition> partitions) {
        List<Double> approximateShapes = new ArrayList<>();

        for (BCCDCladePartition partition : partitions) {
            if (partition.getNumberOfOccurrences() < 2) continue;
            double[] b = partition.getObservedBranchLengthsOld().toArray();
            double approximateShape = GammaDistribution.estimateShapeMLE(b, 1e-6);
            approximateShapes.add(approximateShape);
        }

        Median median = new Median();
        return median.evaluate(approximateShapes.stream().mapToDouble(x -> x).toArray());
    }

    private class BetaOptimizer implements UnivariateFunction {
        BCCD bccd;
        List<BCCDCladePartition> partitions;
        Function<CladePartitionObservation, Double>[] getBetaObservations;
        int[] betaGroup;
        double approximateShape;

        public BetaOptimizer(
                BCCD bccd,
                int[] betaGroup,
                Function<CladePartitionObservation, Double>[] getBetaObservations,
                double approximateShape
        ) {
            this.bccd = bccd;
            this.betaGroup = betaGroup;
            this.getBetaObservations = getBetaObservations;
            this.partitions = bccd.getAllPartitions();
            this.approximateShape = approximateShape;
        }

        @Override
        public double value(double beta) {
            double[] shapes = this.getOptimalShapes(beta);
            double[] scales = this.getOptimalScales(beta, shapes);
            return this.getGroupLogLikelihood(shapes, scales, beta);
        }

        public double[] getOptimalShapes(double beta) {
            double[] shapes = new double[this.betaGroup.length];

            for (int i = 0; i < this.betaGroup.length; i++) {
                int pIdx = this.betaGroup[i];
                BCCDCladePartition partition = partitions.get(pIdx);

                if (partition.getNumberOfOccurrences() == 1) {
                    shapes[i] = this.approximateShape;
                    continue;
                }

                double[] transformedObservations = partition.getObservations().stream().mapToDouble(
                        x -> x.branchLengthOld() * Math.exp(Utils.logOrZero(this.getBetaObservations[pIdx].apply(x)) * beta)
                ).toArray();

                shapes[i] = GammaDistribution.estimateShapeMLE(transformedObservations, 1e-6);
            }

            return shapes;
        }

        public double[] getOptimalScales(double beta, double[] shapes) {
            double[] scales = new double[this.betaGroup.length];

            for (int i = 0; i < this.betaGroup.length; i++) {
                int pIdx = this.betaGroup[i];
                BCCDCladePartition partition = partitions.get(pIdx);

                double[] transformedObservations = partition.getObservations().stream().mapToDouble(
                        x -> x.branchLengthOld() * Math.exp(Utils.logOrZero(this.getBetaObservations[pIdx].apply(x)) * beta)
                ).toArray();

                scales[i] = GammaDistribution.estimateScaleMLE(transformedObservations, shapes[i]);
            }

            return scales;
        }

        public double getGroupLogLikelihood(double[] shapes, double[] scales, double beta) {
            double logLikelihood = 0.0;

            for (int i = 0; i < betaGroup.length; i++) {
                int pIdx = betaGroup[i];
                BCCDCladePartition partition = this.partitions.get(pIdx);

                setDistribution(bccd, pIdx, beta, shapes[i], scales[i]);

                for (CladePartitionObservation observation : partition.getObservations()) {
                    logLikelihood += partition.getDistributionFunc().apply(observation).logDensity(observation.branchLengthOld());
                }
            }

            return logLikelihood;
        }
    }

    private void setDistribution(BCCD bccd, int partitionIdx, double beta, double shape, double scale) {
        Function<CladePartitionObservation, Double> getBetaObservation = this.getBetaObservations[partitionIdx];
        BCCDCladePartition partition = bccd.getAllPartitions().get(partitionIdx);
        partition.setDistributionFunc(
                x -> new GammaDistribution(
                        shape,
                        scale / Math.exp(Utils.logOrZero(getBetaObservation.apply(x)) * beta)
                )
        );
    }
}
