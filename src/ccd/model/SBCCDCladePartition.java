package ccd.model;

import beast.base.evolution.tree.Node;

import java.util.LinkedList;
import java.util.List;
import java.util.function.Function;

public class SBCCDCladePartition extends CladePartition {
    private List<SBCCDCladePartitionObservation> observations = new LinkedList<>();

    public SBCCDCladePartition(SBCCDClade parentClade, SBCCDClade[] childClades) {
        super(parentClade, childClades);
    }

    /* -- BOOK KEEPING - BOOK KEEPING -- */

    public static double getMaxDistanceToLeaf(Node vertex) {
        if (vertex.isLeaf()) return 0.0;

        double firstHeight = getBranchLength(vertex, 0) + getMaxDistanceToLeaf(vertex.getChild(0));
        double secondHeight = getBranchLength(vertex, 1) + getMaxDistanceToLeaf(vertex.getChild(1));

        return Math.max(firstHeight, secondHeight);
    }

    protected static double getBranchLength(Node vertex, int childIndex) {
        double vertexHeight = vertex.getHeight();
        double childHeight = vertex.getChild(childIndex).getHeight();

        if (vertexHeight - childHeight < 0)
            throw new IllegalArgumentException();

        return vertexHeight - childHeight;
    }

    private static SBCCDCladePartitionObservation createObservation(Node vertex) {
        SBCCDCladePartitionObservation observation = new SBCCDCladePartitionObservation(
                getBranchLength(vertex, Utils.getFirstChild(vertex)),
                getBranchLength(vertex, Utils.getSecondChild(vertex)),
                getMaxDistanceToLeaf(vertex),
                getMaxDistanceToLeaf(vertex.getTree().getRoot())
        );
        return observation;
    }

    @Override
    protected void increaseOccurrenceCount(Node vertex) {
        super.increaseOccurrenceCount(vertex);

        SBCCDCladePartitionObservation observation = createObservation(vertex);
        this.observations.add(observation);
    }

    @Override
    protected void decreaseOccurrenceCount(Node vertex) {
        super.decreaseOccurrenceCount(vertex);

        SBCCDCladePartitionObservation observation = createObservation(vertex);
        this.observations.removeIf(x -> x == observation);
    }

    public List<SBCCDCladePartitionObservation> getObservations() {
        return observations;
    }

    /* -- DISTRIBUTION - DISTRIBUTION -- */

    protected double firstBranchAlpha;
    protected double firstBranchBeta;

    protected double secondBranchAlpha;
    protected double secondBranchBeta;

    public void setFirstBranchAlpha(double firstBranchAlpha) {
        this.firstBranchAlpha = firstBranchAlpha;
    }

    public void setFirstBranchBeta(double firstBranchBeta) {
        this.firstBranchBeta = firstBranchBeta;
    }

    public void setSecondBranchAlpha(double secondBranchAlpha) {
        this.secondBranchAlpha = secondBranchAlpha;
    }

    public void setSecondBranchBeta(double secondBranchBeta) {
        this.secondBranchBeta = secondBranchBeta;
    }

    public BranchLengthDistribution getFirstBranchDistribution() {
        return new BetaDistribution(this.firstBranchAlpha, this.firstBranchBeta);
    }

    public BranchLengthDistribution getSecondBranchDistribution() {
        return new BetaDistribution(this.secondBranchAlpha, this.secondBranchBeta);
    }

    public double getCCP(Node vertex) {
        double ccdCCP = super.getCCP();

        double subTreeHeight = this.getMaxDistanceToLeaf(vertex);

        double firstBranchLength = getBranchLength(vertex, Utils.getFirstChild(vertex));
        double secondBranchLength = getBranchLength(vertex, Utils.getSecondChild(vertex));

        BranchLengthDistribution firstBranchLengthDistribution = this.getFirstBranchDistribution();
        BranchLengthDistribution secondBranchLengthDistribution = this.getSecondBranchDistribution();

        double firstBranchProbability;
        if (vertex.getChild(Utils.getFirstChild(vertex)).isLeaf()) {
            firstBranchProbability = 1.0;
        } else {
            firstBranchProbability = firstBranchLengthDistribution.density(firstBranchLength / subTreeHeight) / subTreeHeight;
        }

        double secondBranchProbability;
        if (vertex.getChild(Utils.getSecondChild(vertex)).isLeaf()) {
            secondBranchProbability = 1.0;
        } else {
            secondBranchProbability = secondBranchLengthDistribution.density(secondBranchLength / subTreeHeight) / subTreeHeight;
        }

        return ccdCCP * firstBranchProbability * secondBranchProbability;
    }

    public double getLogCCP(Node vertex) {
        double logCcdCCP = super.getLogCCP();

        double subTreeHeight = this.getMaxDistanceToLeaf(vertex);

        double firstBranchLength = getBranchLength(vertex, Utils.getFirstChild(vertex));
        double secondBranchLength = getBranchLength(vertex, Utils.getSecondChild(vertex));

        BranchLengthDistribution firstBranchLengthDistribution = this.getFirstBranchDistribution();
        BranchLengthDistribution secondBranchLengthDistribution = this.getSecondBranchDistribution();

        double firstLogBranchProbability;
        if (vertex.getChild(Utils.getFirstChild(vertex)).isLeaf()) {
            firstLogBranchProbability = 0.0;
        } else {
            // TODO: do we need a correction factor - log(subTreeHeight) here?
            firstLogBranchProbability = firstBranchLengthDistribution.logDensity(firstBranchLength / subTreeHeight);
        }

        double secondLogBranchProbability;
        if (vertex.getChild(Utils.getSecondChild(vertex)).isLeaf()) {
            secondLogBranchProbability = 0.0;
        } else {
            // TODO: do we need a correction factor - log(subTreeHeight) here?
            secondLogBranchProbability = secondBranchLengthDistribution.logDensity(secondBranchLength / subTreeHeight);
        }

        return logCcdCCP + firstLogBranchProbability + secondLogBranchProbability;
    }

    /* -- Sampling - Sampling -- */

    public double sampleFirstBranchLength(double subTreeHeight) {
        BranchLengthDistribution dist = this.getFirstBranchDistribution();
        double sample = subTreeHeight * dist.sample();
        return sample;
    }

    public double sampleSecondBranchLength(double subTreeHeight) {
        BranchLengthDistribution dist = this.getSecondBranchDistribution();
        double sample = subTreeHeight * dist.sample();
        return sample;
    }

    /* -- MAP Tree - MAP Tree -- */

    public double getMAPFirstBranchLength(double subTreeHeight) {
        BranchLengthDistribution dist = this.getFirstBranchDistribution();
        try {
            return subTreeHeight * dist.mode();
        } catch (NoModeException e) {
            return subTreeHeight * dist.mean();
        }
    }

    public double getMAPSecondBranchLength(double subTreeHeight) {
        BranchLengthDistribution dist = this.getSecondBranchDistribution();
        try {
            return subTreeHeight * dist.mode();
        } catch (NoModeException e) {
            return subTreeHeight * dist.mean();
        }
    }
}
