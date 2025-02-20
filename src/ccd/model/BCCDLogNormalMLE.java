package ccd.model;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.util.*;
import java.util.function.Function;
import java.util.stream.DoubleStream;


public class BCCDLogNormalMLE extends ParameterEstimator<BCCD> {
    Function<CladePartitionObservation, Double>[] getBetaObservations;
    Function<CladePartitionObservation, Double> func;
    BCCDFeatureSelector featureSelector;
    int[][] betaGroups;

    public BCCDLogNormalMLE(Function<CladePartitionObservation, Double> func) {
        this.func = func;
    }

    public BCCDLogNormalMLE(BCCDFeatureSelector featureSelector) {
        this.featureSelector = featureSelector;
    }

    @Override
    public BCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new BCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(BCCD bccd) {
        List<BCCDCladePartition> partitions = bccd.getAllPartitions();

        if (this.func != null) {
            this.getBetaObservations = new Function[partitions.size()];
            for (int i = 0; i < partitions.size(); i++) {
                this.getBetaObservations[i] = func;
            }
        } else {
            this.getBetaObservations = this.featureSelector.getFeatures(partitions);
        }

        this.betaGroups = new int[partitions.size()][];
        for (int i = 0; i < partitions.size(); i++) {
            this.betaGroups[i] = new int[]{i};
        }

        double approximateSigma = this.estimateMedianApproximateSigma(partitions);

        for (int[] betaGroup : this.betaGroups) {
            BrentOptimizer solver = new BrentOptimizer(1e-3, 1e-3);
            BetaOptimizer optimizer = new BetaOptimizer(
                    bccd,
                    betaGroup,
                    this.getBetaObservations,
                    approximateSigma
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

            double[] optimalMus = optimizer.getOptimalMus(optimalBeta);
            double[] optimalSigmas = optimizer.getOptimalSigmas(optimalBeta);

            for (int i = 0; i < betaGroup.length; i++) {
                setDistribution(
                        partitions.get(betaGroup[i]), betaGroup[i], optimalBeta, optimalMus[i], optimalSigmas[i]
                );
            }
        }
    }

    @Override
    public int getNumberOfParameters(BCCD ccd) {
        return 3 * ccd.getNumberOfCladePartitions();
    }

    private double estimateMedianApproximateSigma(List<BCCDCladePartition> partitions) {
        List<Double> approximateSigmas = new ArrayList<>();

        for (BCCDCladePartition partition : partitions) {
            if (partition.getNumberOfOccurrences() < 5) continue;
            double[] b = partition.getObservedLogBranchLengthsOld().toArray();
            double approximateSigma = new Variance().evaluate(b);
            approximateSigmas.add(approximateSigma);
        }

        Median median = new Median();
        return median.evaluate(approximateSigmas.stream().mapToDouble(x -> x).toArray());
    }

    private class BetaOptimizer implements UnivariateFunction {
        BCCD bccd;
        List<BCCDCladePartition> partitions;
        Function<CladePartitionObservation, Double>[] getBetaObservations;
        int[] betaGroup;
        double approximateSigma;

        public BetaOptimizer(
                BCCD bccd,
                int[] betaGroup,
                Function<CladePartitionObservation, Double>[] getBetaObservations,
                double approximateSigma
        ) {
            this.bccd = bccd;
            this.betaGroup = betaGroup;
            this.getBetaObservations = getBetaObservations;
            this.partitions = bccd.getAllPartitions();
            this.approximateSigma = approximateSigma;
        }

        @Override
        public double value(double beta) {
            double[] mus = this.getOptimalMus(beta);
            double[] sigmas = this.getOptimalSigmas(beta);
            return this.getGroupLogLikelihood(mus, sigmas, beta);
        }

        public double[] getOptimalMus(double beta) {
            double[] mus = new double[this.betaGroup.length];

            for (int i = 0; i < this.betaGroup.length; i++) {
                int pIdx = this.betaGroup[i];
                BCCDCladePartition partition = partitions.get(pIdx);

                mus[i] = partition.getObservedLogBranchLengthsOld().average().orElseThrow();

                DoubleStream betaObservations = partition.getObservations().stream().mapToDouble(x -> this.getBetaObservations[pIdx].apply(x));
                mus[i] -= beta * betaObservations.average().orElseThrow();
            }

            return mus;
        }

        public double[] getOptimalSigmas(double beta) {
            double[] sigmas = new double[this.betaGroup.length];

            for (int i = 0; i < this.betaGroup.length; i++) {
                int pIdx = this.betaGroup[i];
                BCCDCladePartition partition = partitions.get(pIdx);

                if (partition.getNumberOfOccurrences() < 5) {
                    sigmas[i] = this.approximateSigma;
                    continue;
                }

                double[] b = partition.getObservedLogBranchLengthsOld().toArray();
                DoubleStream betaObservations = partition.getObservations().stream().mapToDouble(x -> this.getBetaObservations[pIdx].apply(x));
                sigmas[i] = new Variance().evaluate(b) - Math.pow(beta, 2) * new Variance().evaluate(betaObservations.toArray());

                if (sigmas[i] <= 0.0) {
                    sigmas[i] = this.approximateSigma;
                }
            }

            return sigmas;
        }

        public double getGroupLogLikelihood(double[] mus, double[] sigmas, double beta) {
            double logLikelihood = 0.0;

            for (int i = 0; i < betaGroup.length; i++) {
                int pIdx = betaGroup[i];
                BCCDCladePartition partition = this.partitions.get(pIdx);

                setDistribution(partition, pIdx, beta, mus[i], sigmas[i]);

                for (CladePartitionObservation observation : partition.getObservations()) {
                    logLikelihood += partition.getDistributionFunc().apply(observation).logDensity(observation.branchLengthOld());
                }
            }

            return logLikelihood;
        }
    }

    private void setDistribution(BCCDCladePartition partition, int partitionIdx, double beta, double mu, double sigma) {
        Function<CladePartitionObservation, Double> getBetaObservation = this.getBetaObservations[partitionIdx];
        partition.setDistributionFunc(
                x -> new LogNormalDistribution(
                        mu + beta * getBetaObservation.apply(x),
                        Math.sqrt(sigma)
                )
        );
    }
}
