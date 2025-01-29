package ccd.model;

import beast.base.evolution.tree.Tree;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.IllinoisSolver;
import org.apache.commons.math3.exception.NoBracketingException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.DoubleStream;


public class CTMCCCDEstimator extends ParameterEstimator<CTMCCCD> {

    @Override
    public CTMCCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new CTMCCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(CTMCCCD sbccd) {
        this.estimateTimeSinceLastSplitDistribution(sbccd);
        this.estimatePartitionRates(sbccd);
    }

    @Override
    public int getNumberOfParameters(CTMCCCD ccd) {
        return ccd.getNumberOfCladePartitions();
    }

    private void estimateTimeSinceLastSplitDistribution(CTMCCCD sbccd) {
        double[] timesSinceLastSplit = new double[sbccd.getNumberOfBaseTrees()];

        for (int i = 0; i < sbccd.getNumberOfBaseTrees(); i++) {
            Tree tree = sbccd.getBaseTrees().get(i);
            timesSinceLastSplit[i] = Math.log(Utils.getTimeSinceLastSplit(tree.getRoot()));
        }

        BranchLengthDistribution distribution = LogNormalDistribution.estimateMLE(
                timesSinceLastSplit
        );
        sbccd.setTimeSinceLastSplitDistribution(distribution);
    }


    private void estimatePartitionRates(CTMCCCD sbccd) {
        for (CTMCCCDCladePartition partition : sbccd.getAllPartitions()) {
            double numObservations = partition.getNumberOfOccurrences();

            double totalObservedBranchLength = 0.0;
            for (CladePartition parentPartition : partition.getParentClade().getPartitions()) {
                totalObservedBranchLength += ((CTMCCCDCladePartition) parentPartition).getObservations().stream().mapToDouble(CTMCCCDCladePartitionObservation::branchLengthTop).sum();
            }

            double estimatedRate = numObservations / totalObservedBranchLength;
            partition.setRate(estimatedRate);
        }
    }

}
