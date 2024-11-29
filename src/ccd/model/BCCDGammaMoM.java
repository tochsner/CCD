package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;


public class BCCDGammaMoM extends ParameterEstimator<BCCD> {
    List<Function<CladePartitionObservation, Double>> getObservation;
    int numBetas;
    boolean useGlobalBeta;

    public BCCDGammaMoM() {
        this(new ArrayList<>(), true);
    }

    public BCCDGammaMoM(
            Function<CladePartitionObservation, Double> getObservation,
            boolean useGlobalBeta
    ) {
        this(List.of(getObservation), useGlobalBeta);
    }

    public BCCDGammaMoM(
            List<Function<CladePartitionObservation, Double>> getObservation,
            boolean useGlobalBeta
    ) {
        this.getObservation = getObservation;
        this.useGlobalBeta = useGlobalBeta;
        this.numBetas = this.getObservation.size();
    }

    @Override
    public BCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new BCCD(numLeaves, storeBaseTrees, this);
    }

    DoubleStream getLogObservations(int betaIdx, BCCDCladePartition partition) {
        return partition.getObservations().stream().mapToDouble(x -> Utils.logOrZero(this.getObservation.get(betaIdx).apply(x)));
    }

    DoubleStream getObservations(int betaIdx, BCCDCladePartition partition) {
        return partition.getObservations().stream().mapToDouble(x -> this.getObservation.get(betaIdx).apply(x));
    }

    @Override
    public void estimateParameters(BCCD bccd) {
        List<BCCDCladePartition> partitions = bccd.getAllPartitions();

        double[][] betas = getBetas(partitions);
        double[] shapes = getShapes(partitions, betas);
        double[] scales = getScales(partitions, betas, shapes);

        for (int i = 0; i < partitions.size(); i++) {
            int pIdx = i;
            double shape = shapes[i];
            double scale = scales[i];

            BCCDCladePartition partition = partitions.get(i);
            partition.setDistributionFunc(
                    x -> new GammaDistribution(
                            shape,
                            scale * Math.exp(IntStream.range(0, this.numBetas).mapToDouble(j ->  Utils.logOrZero(this.getObservation.get(j).apply(x)) * betas[pIdx][j]).sum())
                    )
            );
        }
    }

    public double[][] getBetas(List<BCCDCladePartition> partitions) {
        double[][] betas = new double[partitions.size()][this.numBetas];

        for (int i = 0; i < this.numBetas; i++) {
            double beta = 0.0;
            int numObservations = 0;

            for (int j = 0; j < partitions.size(); j++) {
                BCCDCladePartition partition = partitions.get(j);

                double[] b1 = partition.getObservedLogBranchLengthsOld().toArray();
                double[] b2 = this.getLogObservations(i, partition).toArray();

                if (b1.length < 2) {
                    continue;
                }

                double partitionBeta = new Covariance().covariance(b1, b2) / new Variance().evaluate(b2);

                if (Double.isNaN(partitionBeta)) {
                    continue;
                }

                if (this.useGlobalBeta) {
                    beta += b1.length * partitionBeta;
                    numObservations += b1.length;
                } else {
                    betas[j][i] = partitionBeta;
                }
            }

            if (this.useGlobalBeta) {
                for (int j = 0; j < partitions.size(); j++) {
                    betas[j][i] = beta / numObservations;
                }
            }
        }

        return betas;
    }

    public double[] getShapes(List<BCCDCladePartition> partitions, double[][] betas) {
        double[] shapes = new double[partitions.size()];
        List<Double> allShapes = new LinkedList<>();

        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Integer> idxWithNegativeEstimate = new ArrayList<>();

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);
            int pIdx = i;

            double[] b = partition.getObservedLogBranchLengthsOld().toArray();
            if (b.length < 2) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            double trigammaShape = new Variance().evaluate(b);
            for (int j = 0; j < this.numBetas; j++) {
                trigammaShape -= Math.pow(betas[pIdx][j], 2) * new Variance().evaluate(this.getLogObservations(j, partition).toArray());
            }

            if (trigammaShape <= 0.0) {
                idxWithNegativeEstimate.add(i);
                continue;
            }
            double shape = InverseTrigamma.value(trigammaShape);

            if (shape < 0) {
                idxWithNegativeEstimate.add(i);
                continue;
            }

            shapes[i] = shape;
            allShapes.add(shape);
        }

        double meanShape = allShapes.stream().reduce((a, b) -> a + b).get() / partitions.size();
        for (int i : idxWithoutEnoughData) {
            shapes[i] = meanShape;
        }

        double lowShape = new Percentile().evaluate(allShapes.stream().mapToDouble(x -> x).toArray(), 10);
        for (int i : idxWithNegativeEstimate) {
            shapes[i] = lowShape;
        }

        return shapes;
    }

    public double[] getScales(List<BCCDCladePartition> partitions, double[][] betas, double[] shapes) {
        double[] scales = new double[partitions.size()];

        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Double> allScales = new LinkedList<>();

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double scale = partition.getObservedBranchLengthsOld().average().orElseThrow() / shapes[i];

            for (int j = 0; j < this.numBetas; j++) {
                final int pIdx = i;
                final int pBeta = j;
                scale /= this.getObservations(j, partition).map(x -> Math.exp(betas[pIdx][pBeta] * Utils.logOrZero(x))).average().orElseThrow();
            }

            if (scale == 0.0) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            scales[i] = scale;
            allScales.add(scale);
        }

        double meanScale = allScales.stream().reduce((a, b) -> a + b).get() / partitions.size();
        for (int i : idxWithoutEnoughData) {
            scales[i] = meanScale;
        }

        return scales;
    }
}
