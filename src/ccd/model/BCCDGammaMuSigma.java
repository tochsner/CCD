package ccd.model;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class BCCDGammaMuSigma extends ParameterEstimator<BCCD> {
    @Override
    public BCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new BCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(BCCD bccd) {
        List<BCCDCladePartition> partitions = bccd.getAllPartitions();
        double[] scales = getScales(partitions);
        double[] shapes = getShapes(partitions, scales);

        for (int i = 0; i < partitions.size(); i++) {
            double shape = shapes[i];
            double scale = scales[i];

            BCCDCladePartition partition = partitions.get(i);
            partition.setDistributionFunc(
                    x -> new GammaDistribution(shape, scale)
            );
        }
    }

    public double[] getScales(List<BCCDCladePartition> partitions) {
        double[] scales = new double[partitions.size()];

        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Double> allScales = new LinkedList<>();

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            if (partition.getNumberOfOccurrences() < 2) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            double b = partition.getObservedBranchLengthsOld().average().orElseThrow();
            double logB = partition.getObservedLogBranchLengthsOld().average().orElseThrow();
            double bLogB = partition.getObservedBranchLengthsOld().map(x -> x * Math.log(x)).average().orElseThrow();

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

    public double[] getShapes(List<BCCDCladePartition> partitions, double[] scales) {
        double[] shapes = new double[partitions.size()];

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double b = partition.getObservedBranchLengthsOld().average().orElseThrow();

            shapes[i] = b / scales[i];
        }

        return shapes;
    }
}
