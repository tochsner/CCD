package ccd.model;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;


public class BCCDGammaLinear extends BCCDParameterEstimator {
    List<Function<BCCDCladePartition, DoubleStream>> getObservations;
    List<Function<CladePartitionObservation, Double>> getObservation;
    int numBetas;
    boolean useGlobalBeta;

    public BCCDGammaLinear(
            List<Function<BCCDCladePartition, DoubleStream>> getObservations,
            List<Function<CladePartitionObservation, Double>> getObservation,
            boolean useGlobalBeta
    ) {
        this.getObservation = getObservation;
        this.getObservations = getObservations;
        if (this.getObservation.size() != this.getObservations.size())
            throw new IllegalArgumentException("Function array lengths must match.");

        this.useGlobalBeta = useGlobalBeta;
        this.numBetas = this.getObservation.size();
    }

    @Override
    public void estimateParameters(List<BCCDCladePartition> partitions) {
        double[][] betas = getBetas(partitions);
        double[] scales = getScales(partitions, betas);
        double[] shapes = getShapes(partitions, scales, betas);

        for (int i = 0; i < partitions.size(); i++) {
            double shape = shapes[i];
            double scale = scales[i];

            BCCDCladePartition partition = partitions.get(i);
            partition.setDistributionFunc(
                    x -> new GammaDistribution(shape, scale)
            );
        }
    }

    public double[][] getBetas(List<BCCDCladePartition> partitions) {
        return new double[partitions.size()][this.numBetas];
    }

    public double[] getScales(List<BCCDCladePartition> partitions, double[][] betas) {
        double[] scales = new double[partitions.size()];

        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Double> allScales = new LinkedList<>();

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            if (partition.getNumberOfOccurrences() < 2) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            double b = partition.getObservedLogBranchLengthsOld().average().orElseThrow();
            double logB = partition.getObservedLogBranchLengthsOld().map(x -> Math.log(x)).average().orElseThrow();
            double bLogB = partition.getObservedLogBranchLengthsOld().map(x -> x * Math.log(x)).average().orElseThrow();

            double scale = bLogB - b * logB;

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

    public double[] getShapes(List<BCCDCladePartition> partitions, double[] scales, double[][] betas) {
        double[] shapes = new double[partitions.size()];

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double b = partition.getObservedLogBranchLengthsOld().average().orElseThrow();

            shapes[i] = b / scales[i];
        }

        return shapes;
    }
}
