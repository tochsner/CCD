package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.function.Function;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;


public class BCCDLinear extends BCCDParameterEstimator {
    List<Function<BCCDCladePartition, DoubleStream>> getObservations;
    List<Function<CladePartitionObservation, Double>> getObservation;

    public BCCDLinear(
            List<Function<BCCDCladePartition, DoubleStream>> getObservations,
            List<Function<CladePartitionObservation, Double>> getObservation
    ) {
        this.getObservation = getObservation;
        this.getObservations = getObservations;

        if (this.getObservation.size() != this.getObservations.size())
            throw new IllegalArgumentException("Function array lengths must match.");
    }

    @Override
    public void estimateParameters(List<BCCDCladePartition> partitions) {
        double[] betas = getBetas(partitions);
        double[] mus = getMus(partitions, betas);
        double[] sigmas = getSigmas(partitions, betas);

        for (int i = 0; i < partitions.size(); i++) {
            double mu = mus[i];
            double sigma = sigmas[i];

            BCCDCladePartition partition = partitions.get(i);
            partition.setLogMeanFunc(x -> mu + IntStream.range(0, betas.length).mapToDouble(j -> betas[j] * this.getObservation.get(j).apply(x)).sum());
            partition.setLogVarianceFunc(x -> sigma);
        }
    }

    public double[] getBetas(List<BCCDCladePartition> partitions) {
        double[] betas = new double[this.getObservation.size()];

        for (int i = 0; i < this.getObservation.size(); i++) {
            int numObservations = 0;

            for (int j = 0; j < partitions.size(); j++) {
                BCCDCladePartition partition = partitions.get(j);

                double[] b1 = partition.getObservedLogBranchLengthsOld().toArray();
                double[] b2 = this.getObservations.get(i).apply(partition).toArray();

                if (b1.length < 2) {
                    continue;
                }

                double partitionBeta = new Covariance().covariance(b1, b2) / new Variance().evaluate(b2);

                if (Double.isNaN(partitionBeta)) {
                    continue;
                }

                betas[i] += b1.length * partitionBeta;
                numObservations += b1.length;
            }

            betas[i] /= numObservations;
        }

        return betas;
    }

    public double[] getMus(List<BCCDCladePartition> partitions, double[] betas) {
        double[] mus = new double[partitions.size()];

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            mus[i] = partition.getObservedLogBranchLengthsOld().average().orElseThrow();

            for (int j = 0; j < this.getObservation.size(); j++) {
                mus[i] -= betas[j] * this.getObservations.get(j).apply(partition).average().orElseThrow();
            }
        }

        return mus;
    }

    public double[] getSigmas(List<BCCDCladePartition> partitions, double[] betas) {
        double[] sigmas = new double[partitions.size()];

        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Integer> idxWithNegativeEstimate = new ArrayList<>();

        List<Double> allSigmas = new LinkedList<>();

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double[] b = partition.getObservedLogBranchLengthsOld().toArray();
            if (b.length < 2) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            double sigma = new Variance().evaluate(b);
            for (int j = 0; j < this.getObservation.size(); j++) {
                 sigma -= Math.pow(betas[j], 2) * new Variance().evaluate(this.getObservations.get(j).apply(partition).toArray());
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
}
