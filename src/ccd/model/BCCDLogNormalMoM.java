package ccd.model;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.function.Function;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;


public class BCCDLogNormalMoM extends ParameterEstimator<BCCD> {
    List<Function<CladePartitionObservation, Double>> getObservation;
    int numBetas;
    boolean useGlobalBeta;

    public BCCDLogNormalMoM() {
        this(new ArrayList<>(), true);
    }

    public BCCDLogNormalMoM(
            Function<CladePartitionObservation, Double> getObservation,
            boolean useGlobalBeta
    ) {
        this(List.of(getObservation), useGlobalBeta);
    }

    public BCCDLogNormalMoM(
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

    DoubleStream getObservations(int betaIdx, BCCDCladePartition partition) {
        return partition.getObservations().stream().mapToDouble(x -> this.getObservation.get(betaIdx).apply(x));
    }

    @Override
    public void estimateParameters(BCCD bccd) {
        List<BCCDCladePartition> partitions = bccd.getAllPartitions();

        double[][] betas = this.getBetas(partitions);
        double[] mus = this.getMus(partitions, betas);
        double[] sigmas = this.getSigmas(partitions, betas);

        double[] modes = this.getModes(partitions, betas, mus, sigmas);

        for (int i = 0; i < partitions.size(); i++) {
            double mu = mus[i];
            double sigma = sigmas[i];
            int pIdx = i;

            BCCDCladePartition partition = partitions.get(i);
            partition.setDistributionFunc(
                    x -> new LogNormalDistribution(
                            mu + IntStream.range(0, this.numBetas).mapToDouble(j -> betas[pIdx][j] * this.getObservation.get(j).apply(x)).sum(),
                            Math.sqrt(sigma)
                    )
            );
        }
    }

    @Override
    public int getNumberOfParameters(BCCD ccd) {
        return 2 * ccd.getNumberOfCladePartitions() + numBetas * (this.useGlobalBeta ? 1 : ccd.getNumberOfCladePartitions());
    }

    public double[][] getBetas(List<BCCDCladePartition> partitions) {
        double[][] betas = new double[partitions.size()][this.numBetas];

        for (int i = 0; i < this.numBetas; i++) {
            double beta = 0.0;
            int numObservations = 0;

            for (int j = 0; j < partitions.size(); j++) {
                BCCDCladePartition partition = partitions.get(j);

                double[] b1 = partition.getObservedLogBranchLengthsOld().toArray();
                double[] b2 = this.getObservations(i, partition).toArray();

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

    public double[] getMus(List<BCCDCladePartition> partitions, double[][] betas) {
        double[] mus = new double[partitions.size()];

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);
            int pIdx = i;

            mus[i] = partition.getObservedLogBranchLengthsOld().average().orElseThrow();

            for (int j = 0; j < this.numBetas; j++) {
                mus[i] -= betas[pIdx][j] * this.getObservations(j, partition).average().orElseThrow();
            }
        }

        return mus;
    }

    public double[] getSigmas(List<BCCDCladePartition> partitions, double[][] betas) {
        double[] sigmas = new double[partitions.size()];

        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Integer> idxWithNegativeEstimate = new ArrayList<>();

        List<Double> allSigmas = new LinkedList<>();

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);
            int pIdx = i;

            double[] b = partition.getObservedLogBranchLengthsOld().toArray();
            if (b.length < 2) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            double sigma = new Variance().evaluate(b);
            for (int j = 0; j < this.numBetas; j++) {
                sigma -= Math.pow(betas[pIdx][j], 2) * new Variance().evaluate(this.getObservations(j, partition).toArray());
            }

            if (sigma <= 0.0) {
                idxWithNegativeEstimate.add(i);
                continue;
            }

            sigmas[i] = sigma;
            allSigmas.add(sigma);
        }

        double meanSigma = allSigmas.stream().reduce((a, b) -> a + b).get() / partitions.size();
        for (int i : idxWithoutEnoughData) {
            sigmas[i] = meanSigma;
        }

        double lowSigma = new Percentile().evaluate(allSigmas.stream().mapToDouble(x -> x).toArray(), 10);
        for (int i : idxWithNegativeEstimate) {
            sigmas[i] = lowSigma;
        }

        return sigmas;
    }

    public double[] getModes(List<BCCDCladePartition> partitions, double[][] betas, double[] mus, double[] sigmas) {
        // build up linear system

        int n = partitions.size();
        RealMatrix A = new BlockRealMatrix(n, n);
        RealVector b = new ArrayRealVector(n);

        int rootPartition = IntStream.range(0, partitions.size()).filter(i -> partitions.get(i).getParentClade().isRoot()).findFirst().orElseThrow();
        this.buildCoefficientMatrix(rootPartition, -1, A, b, partitions, betas, mus, sigmas);

        for (BCCDCladePartition partition : partitions) {
            // pass
        }

        return null;
    }

    private void buildCoefficientMatrix(
            int partitionIdx,
            int parentPartitionIdx,
            RealMatrix A,
            RealVector b,
            List<BCCDCladePartition> partitions,
            double[][] betas,
            double[] mus,
            double[] sigmas
    ) {
        if (parentPartitionIdx != -1) {
            // take over parent coefficients

            for (int i = 0; i < partitions.size(); i++) {
                A.addToEntry(partitionIdx, i, A.getEntry(parentPartitionIdx, i) * betas[0][parentPartitionIdx]);
            }

            b.addToEntry(partitionIdx, b.getEntry(parentPartitionIdx) * betas[0][parentPartitionIdx]);
        }

        b.addToEntry(partitionIdx, -mus[partitionIdx] / sigmas[partitionIdx]);

        BCCDCladePartition partition = partitions.get(partitionIdx);

        // pass
    }
}
