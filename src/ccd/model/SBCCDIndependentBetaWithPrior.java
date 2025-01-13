package ccd.model;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.util.LinkedList;
import java.util.List;
import java.util.stream.DoubleStream;


public class SBCCDIndependentBetaWithPrior extends ParameterEstimator<SBCCD> {

    @Override
    public SBCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new SBCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(SBCCD sbccd) {
        this.estimateHeightDistribution(sbccd);
        this.estimatePartitionParameters(sbccd);
    }

    @Override
    public int getNumberOfParameters(SBCCD ccd) {
        int heightParameters = 2;

        int fractionParameters = 0;
        for (SBCCDCladePartition partition : ccd.getAllPartitions()) {
            if (!partition.getChildClades()[0].isLeaf()) fractionParameters++;
            if (!partition.getChildClades()[1].isLeaf()) fractionParameters++;
        }

        return heightParameters + fractionParameters;
    }

    private void estimateHeightDistribution(SBCCD sbccd) {
        Clade rootClade = sbccd.getRootClade();

        DoubleStream observedHeights = DoubleStream.empty();

        for (SBCCDCladePartition partition : sbccd.getAllPartitions()) {
            if (partition.getParentClade() == rootClade) {
                observedHeights = DoubleStream.concat(
                        observedHeights,
                        partition.getObservations().stream().mapToDouble(x -> Utils.logOrZero(x.subTreeHeight()))
                );
            }
        }

        double[] observedLogHeightsArray = observedHeights.toArray();
        LogNormalDistribution heightDistribution = LogNormalDistribution.estimateMLE(observedLogHeightsArray);
        sbccd.setHeightDistribution(heightDistribution);
    }

    private void estimatePartitionParameters(SBCCD sbccd) {
        List<SBCCDCladePartition> partitions = sbccd.getAllPartitions();

        List<Integer> idxWithoutEnoughData = new LinkedList<>();
        List<Double> allAlphas = new LinkedList<>();
        List<Double> allBetas = new LinkedList<>();

        for (int i = 0; i < partitions.size(); i++) {
            SBCCDCladePartition partition = partitions.get(i);

            if (!Utils.getFirstClade(partition).isLeaf()) {
                double[] firstBranchRatios = partition.getObservations().stream().mapToDouble(x -> x.branchLengthLeft() / x.subTreeHeight()).toArray();

                try {
                    double[] params = BetaDistribution.estimateMAP(firstBranchRatios, 1e-2);

                    partition.setFirstBranchAlpha(params[0]);
                    partition.setFirstBranchBeta(params[1]);

                    allAlphas.add(params[0]);
                    allBetas.add(params[1]);
                } catch (Exception e) {
                }
            }

            if (!Utils.getFirstClade(partition).isLeaf()) {
                double[] secondBranchRatios = partition.getObservations().stream().mapToDouble(x -> x.branchLengthRight() / x.subTreeHeight()).toArray();

                try {
                    double[] params = BetaDistribution.estimateMAP(secondBranchRatios, 1e-6);

                    partition.setSecondBranchAlpha(params[0]);
                    partition.setSecondBranchBeta(params[1]);

                    allAlphas.add(params[0]);
                    allBetas.add(params[1]);
                } catch (Exception e) {
                }
            }
        }

        Median median = new Median();
        double medianAlpha = median.evaluate(allAlphas.stream().mapToDouble(x -> x).toArray());
        double medianBeta = median.evaluate(allBetas.stream().mapToDouble(x -> x).toArray());

        for (int i : idxWithoutEnoughData) {
            SBCCDCladePartition partition = partitions.get(i);

            if (!Utils.getFirstClade(partition).isLeaf()) {
                partition.setFirstBranchAlpha(medianAlpha);
                partition.setFirstBranchBeta(medianBeta);
            }

            if (!Utils.getSecondClade(partition).isLeaf()) {
                partition.setSecondBranchAlpha(medianAlpha);
                partition.setSecondBranchBeta(medianBeta);
            }
        }
    }

}
