package ccd.model;

import beast.base.evolution.tree.Node;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.ConstantRealDistribution;

import java.util.HashMap;
import java.util.Map;

public class BCCDCladePartition extends CladePartition {
    private Map<Double, Integer> minBranchLengthOccurrences = new HashMap<>();

    public BCCDCladePartition(BCCDClade parentClade, BCCDClade[] childClades) {
        super(parentClade, childClades);
    }

    /* -- BOOK KEEPING - BOOK KEEPING -- */

    protected static double getMinLogBranchLength(Node vertex) {
        if (vertex.isLeaf()) {
            return 0.0;
        }

        double vertexHeight = vertex.getHeight();

        double childHeight1 = vertex.getChild(0).getHeight();
        double childHeight2 = vertex.getChild(1).getHeight();

        double branchLength1 = vertexHeight - childHeight1;
        double branchLength2 = vertexHeight - childHeight2;

        return Math.log(Math.min(branchLength1, branchLength2));
    }

    @Override
    protected void increaseOccurrenceCount(Node vertex) {
        super.increaseOccurrenceCount(vertex);

        double branchLength = this.getMinLogBranchLength(vertex);
        this.minBranchLengthOccurrences.merge(branchLength, 1, Integer::sum);
    }

    @Override
    protected void decreaseOccurrenceCount(Node vertex) {
        super.decreaseOccurrenceCount(vertex);

        double branchLength = this.getMinLogBranchLength(vertex);
        this.minBranchLengthOccurrences.merge(branchLength, -1, Integer::sum);
    }

    /* -- CPP - CPP -- */
    @Override
    public double getCCP(Node vertex) {
        double ccdCCP = super.getCCP();

        double minBranchLength = getMinLogBranchLength(vertex);

        AbstractRealDistribution branchLengthDistribution = this.getBranchLengthDistribution(vertex);
        double branchProbability = branchLengthDistribution.density(Math.exp(minBranchLength));

        return ccdCCP * branchProbability;
    }

    /* -- Sampling - Sampling -- */

    protected double getLogMean(Node vertex) {
        double mean = 0;
        int totalSamples = 0;

        for (Map.Entry<Double, Integer> branchLength : this.minBranchLengthOccurrences.entrySet()) {
            mean += branchLength.getValue() * branchLength.getKey();
            totalSamples++;
        }

        mean /= totalSamples;
        return mean;
    }

    protected double getLogVariance(Node vertex, double logSampleMean) {
        double variance = 0;
        int totalSamples = 0;

        for (Map.Entry<Double, Integer> branchLength : this.minBranchLengthOccurrences.entrySet()) {
            variance += branchLength.getValue() * Math.pow(branchLength.getKey() - logSampleMean, 2);
            totalSamples++;
        }

        variance /= totalSamples;

        return variance;
    }

    protected AbstractRealDistribution getBranchLengthDistribution(Node vertex) {
        double scale = this.getLogMean(vertex);
        double shape = this.getLogVariance(vertex, scale);

        if (shape == 0.0) {
            return new ConstantRealDistribution(this.minBranchLengthOccurrences.entrySet().iterator().next().getKey());
        }

        LogNormalDistribution dist = new LogNormalDistribution(scale, shape);
        return dist;
    }

    public double sampleMinBranchLength(Node vertex) {
        AbstractRealDistribution dist = this.getBranchLengthDistribution(vertex);
        return dist.sample();
    }
}
