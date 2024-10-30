package ccd.model;

import beast.base.evolution.tree.Node;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.ConstantRealDistribution;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.stream.DoubleStream;

public class BCCDCladePartition extends CladePartition {
    private List<CladePartitionObservation> observations = new LinkedList<>();

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

    protected static double getMinLogBranchLengthDown(Node vertex) {
        Node firstChild = vertex.getChild(0);
        Node secondChild = vertex.getChild(1);

        double firstChildHeight = firstChild.getHeight();
        double secondChildHeight = secondChild.getHeight();

        if (firstChildHeight < secondChildHeight) {
            return getMinLogBranchLength(secondChild);
        } else {
            return getMinLogBranchLength(firstChild);
        }
    }

    @Override
    protected void increaseOccurrenceCount(Node vertex) {
        super.increaseOccurrenceCount(vertex);

        CladePartitionObservation observation = new CladePartitionObservation(
                getMinLogBranchLength(vertex),
                getMinLogBranchLengthDown(vertex)
        );

        this.observations.add(observation);
    }

    @Override
    protected void decreaseOccurrenceCount(Node vertex) {
        super.decreaseOccurrenceCount(vertex);

        CladePartitionObservation observation = new CladePartitionObservation(
                getMinLogBranchLength(vertex),
                getMinLogBranchLengthDown(vertex)
        );

        this.observations.removeIf(x -> x == observation);
    }

    public List<CladePartitionObservation> getObservations() {
        return observations;
    }

    public DoubleStream getObservedMinLogBranchLength() {
        return this.getObservations().stream().mapToDouble(x -> x.logMinBranchLength());
    }

    public DoubleStream getObservedMinLogBranchLengthsDown() {
        return this.getObservations().stream().mapToDouble(x -> x.logMinBranchLengthDown());
    }

    /* -- DISTRIBUTION - DISTRIBUTION -- */

    protected double logMean;
    protected double logVariance;

    public double getLogMean() {
        return logMean;
    }

    public void setLogMean(double logMean) {
        this.logMean = logMean;
    }

    public double getLogVariance() {
        return logVariance;
    }

    public void setLogVariance(double logVariance) {
        this.logVariance = logVariance;
    }

    @Override
    public double getCCP(Node vertex) {
        double ccdCCP = super.getCCP();

        double minBranchLength = getMinLogBranchLength(vertex);

        AbstractRealDistribution branchLengthDistribution = this.getBranchLengthDistribution(vertex);
        double branchProbability = branchLengthDistribution.density(Math.exp(minBranchLength));

        return ccdCCP * branchProbability;
    }

    /* -- Sampling - Sampling -- */

    protected AbstractRealDistribution getBranchLengthDistribution(Node vertex) {
        double scale = this.getLogMean();
        double shape = this.getLogVariance();

        if (shape == 0.0) {
            return new ConstantRealDistribution(scale);
        } else {
            return new LogNormalDistribution(scale, shape);
        }
    }

    public double sampleMinBranchLength(Node vertex) {
        AbstractRealDistribution dist = this.getBranchLengthDistribution(vertex);
        return dist.sample();
    }
}
