package ccd.model;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.util.LinkedList;
import java.util.List;
import java.util.stream.DoubleStream;


public class SBCCDIndependent extends ParameterEstimator<SBCCD> {

    @Override
    public SBCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new SBCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(SBCCD sbccd) {
        this.estimateHeightDistribution(sbccd);
        this.estimatePartitionParameters(sbccd);
    }

    private void estimateHeightDistribution(SBCCD sbccd) {
        Clade rootClade = sbccd.getRootClade();

        DoubleStream observedLogHeights = DoubleStream.empty();

        for (SBCCDCladePartition partition : sbccd.getAllPartitions()) {
            if (partition.getParentClade() == rootClade) {
                observedLogHeights = DoubleStream.concat(
                        observedLogHeights,
                        partition.getObservations().stream().mapToDouble(x -> Math.log(x.subTreeHeight()))
                );
            }
        }

        double[] observedLogHeightsArray = observedLogHeights.toArray();
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
            if (partition.getNumberOfOccurrences() < 2) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            if (partition.getChildClades()[0].isLeaf()) {
                partition.setFirstBranchDistribution(new UniformDistribution());
            } else {
                double[] firstBranchRatios = partition.getObservations().stream().mapToDouble(x -> x.branchLengthLeft() / x.subTreeHeight()).toArray();

                double alpha = BetaDistribution.estimateAlpha(firstBranchRatios);
                double beta = BetaDistribution.estimateBeta(firstBranchRatios);
                BranchLengthDistribution firstBranchRatioDist = new BetaDistribution(alpha, beta);

                partition.setFirstBranchDistribution(firstBranchRatioDist);

                allAlphas.add(alpha);
                allBetas.add(beta);
            }

            if (partition.getChildClades()[1].isLeaf()) {
                partition.setSecondBranchDistribution(new UniformDistribution());
            } else {
                double[] secondBranchRatios = partition.getObservations().stream().mapToDouble(x -> x.branchLengthRight() / x.subTreeHeight()).toArray();

                double alpha = BetaDistribution.estimateAlpha(secondBranchRatios);
                double beta = BetaDistribution.estimateBeta(secondBranchRatios);
                BranchLengthDistribution secondBranchRatioDist = new BetaDistribution(alpha, beta);

                partition.setSecondBranchDistribution(secondBranchRatioDist);

                allAlphas.add(alpha);
                allBetas.add(beta);
            }
        }

        Median median = new Median();
        double medianAlpha = median.evaluate(allAlphas.stream().mapToDouble(x -> x).toArray());
        double medianBeta = median.evaluate(allBetas.stream().mapToDouble(x -> x).toArray());
        BetaDistribution medianDistribution = new BetaDistribution(medianAlpha, medianBeta);

        for (int i : idxWithoutEnoughData) {
            partitions.get(i).setFirstBranchDistribution(medianDistribution);
            partitions.get(i).setSecondBranchDistribution(medianDistribution);
        }
    }

}
