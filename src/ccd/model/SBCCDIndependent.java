package ccd.model;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

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
        double logMean = new Mean().evaluate(observedLogHeightsArray);
        double logVariance = new Variance().evaluate(observedLogHeightsArray);

        sbccd.setHeightDistribution(
                new LogNormalDistribution(logMean, Math.sqrt(logVariance))
        );
    }

    private void estimatePartitionParameters(SBCCD sbccd) {
        List<SBCCDCladePartition> partitions = sbccd.getAllPartitions();

        for (SBCCDCladePartition partition : partitions) {
            if (partition.getChildClades()[0].isLeaf()) {
                partition.setFirstBranchDistribution(new UniformDistribution());
            } else {
                double[] firstBranchRatios = partition.getObservations().stream().mapToDouble(x -> x.branchLengthLeft() / x.subTreeHeight()).toArray();
                BranchLengthDistribution firstBranchRationDist = BetaDistribution.estimate(firstBranchRatios);
                partition.setFirstBranchDistribution(firstBranchRationDist);
            }

            if (partition.getChildClades()[1].isLeaf()) {
                partition.setSecondBranchDistribution(new UniformDistribution());
            } else {
                double[] secondBranchRatios = partition.getObservations().stream().mapToDouble(x -> x.branchLengthRight() / x.subTreeHeight()).toArray();
                BranchLengthDistribution secondBranchRationDist = BetaDistribution.estimate(secondBranchRatios);
                partition.setSecondBranchDistribution(secondBranchRationDist);
            }
        }
    }

}
