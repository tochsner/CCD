package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.inference.distribution.LogNormalDistributionModel;
import ccd.algorithms.BitSetUtil;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class BCCDCladePartition extends CladePartition {
    private Map<Double, Integer> minBranchLengthOccurrences = new HashMap<>();

    public BCCDCladePartition(BCCDClade parentClade, BCCDClade[] childClades) {
        super(parentClade, childClades);
    }

    /* -- BOOK KEEPING - BOOK KEEPING -- */

    protected static double getMinBranchLength(Node vertex) {
        double vertexHeight = vertex.getHeight();

        double childHeight1 = vertex.getChild(0).getHeight();
        double childHeight2 = vertex.getChild(1).getHeight();

        double branchLength1 = vertexHeight - childHeight1;
        double branchLength2 = vertexHeight - childHeight2;

        return Math.min(branchLength1, branchLength2);
    }

    @Override
    protected void increaseOccurrenceCount(Node vertex) {
        super.increaseOccurrenceCount(vertex);

        double branchLength = this.getMinBranchLength(vertex);
        this.minBranchLengthOccurrences.merge(branchLength, 1, Integer::sum);
    }

    @Override
    protected void decreaseOccurrenceCount(Node vertex) {
        super.decreaseOccurrenceCount(vertex);

        double branchLength = this.getMinBranchLength(vertex);
        this.minBranchLengthOccurrences.merge(branchLength, -1, Integer::sum);
    }

    /* -- Sampling - Sampling -- */

    protected AbstractRealDistribution getBranchLengthDistribution() {
        double scale = 0;
        int totalSamples = 0;
        for (Map.Entry<Double, Integer> branchLength : this.minBranchLengthOccurrences.entrySet()) {
            scale += branchLength.getValue() * Math.log(branchLength.getKey());
            totalSamples++;
        }
        scale /= totalSamples;

        double shape = 0;
        for (Map.Entry<Double, Integer> branchLength : this.minBranchLengthOccurrences.entrySet()) {
            shape += branchLength.getValue() * Math.pow(Math.log(branchLength.getKey()) - scale, 2);
        }
        shape /= totalSamples;

        LogNormalDistribution dist = new LogNormalDistribution(scale, shape);
        return dist;
    }

    public double sampleMinBranchLength() {
        AbstractRealDistribution dist = this.getBranchLengthDistribution();
        return dist.sample();
    }
}
