package ccd.model;

import beast.base.evolution.tree.Node;
import ccd.model.CladePartition;
import ccd.model.CladePartitionObservation;

import java.util.LinkedList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.DoubleStream;

public class BCCDCladePartition extends CladePartition {
    private List<CladePartitionObservation> observations = new LinkedList<>();

    public BCCDCladePartition(BCCDClade parentClade, BCCDClade[] childClades) {
        super(parentClade, childClades);
    }

    /* -- BOOK KEEPING - BOOK KEEPING -- */

    protected static int getOldChildIndex(Node vertex) {
        double childHeight1 = vertex.getChild(0).getHeight();
        double childHeight2 = vertex.getChild(1).getHeight();

        if (childHeight2 < childHeight1) {
            return 0;
        } else {
            return 1;
        }
    }

    protected static int getSmallChildIndex(Node vertex) {
        double childSize1 = vertex.getChild(0).getLeafNodeCount();
        double childSize2 = vertex.getChild(1).getLeafNodeCount();

        if (childSize1 < childSize2) {
            return 0;
        } else {
            return 1;
        }
    }

    protected static double getBranchLength(Node vertex, int childIndex) {
        double vertexHeight = vertex.getHeight();
        double childHeight = vertex.getChild(childIndex).getHeight();

//        if (vertexHeight - childHeight < 0)
//            throw new IllegalArgumentException();

        return vertexHeight - childHeight;
    }

    protected static double getBranchLengthOld(Node vertex) {
        if (vertex.isLeaf()) return 0.0;
        int oldChild = getOldChildIndex(vertex);
        return getBranchLength(vertex, oldChild);
    }

    public static double getMaxDistanceToLeaf(Node vertex) {
        if (vertex.isLeaf()) return 0.0;

        double firstHeight = getBranchLength(vertex, 0) + getMaxDistanceToLeaf(vertex.getChild(0));
        double secondHeight = getBranchLength(vertex, 1) + getMaxDistanceToLeaf(vertex.getChild(1));

        return Math.max(firstHeight, secondHeight);
    }

    protected static double getHeightOld(Node vertex) {
        if (vertex.isLeaf()) return 0.0;
        int oldChild = getOldChildIndex(vertex);
        return getMaxDistanceToLeaf(vertex.getChild(oldChild));
    }

    protected static double getBranchLengthYoung(Node vertex) {
        if (vertex.isLeaf()) return 0.0;
        int youngChild = 1 - getOldChildIndex(vertex);
        return getBranchLength(vertex, youngChild);
    }

    protected static double getBranchLengthOldOld(Node vertex) {
        if (vertex.getLeft().isLeaf() || vertex.getRight().isLeaf()) return 0.0;
        Node oldChild = vertex.getChild(getOldChildIndex(vertex));
        return getBranchLengthOld(oldChild);
    }

    protected static double getBranchLengthOldYoung(Node vertex) {
        if (vertex.getLeft().isLeaf() || vertex.getRight().isLeaf()) return 0.0;
        Node oldChild = vertex.getChild(getOldChildIndex(vertex));
        return getBranchLengthYoung(oldChild);
    }

    protected static double getBranchLengthOldSmall(Node vertex) {
        if (vertex.getLeft().isLeaf() || vertex.getRight().isLeaf()) return 0.0;
        Node oldChild = vertex.getChild(getOldChildIndex(vertex));
        if (oldChild.isLeaf()) return 0.0;
        int smallChild = getSmallChildIndex(oldChild);
        return getBranchLength(oldChild, smallChild);
    }

    protected static double getBranchLengthOldBig(Node vertex) {
        if (vertex.getLeft().isLeaf() || vertex.getRight().isLeaf()) return 0.0;
        Node oldChild = vertex.getChild(getOldChildIndex(vertex));
        if (oldChild.isLeaf()) return 0.0;
        int bigChild = 1 - getSmallChildIndex(oldChild);
        return getBranchLength(oldChild, bigChild);
    }

    private static CladePartitionObservation createObservation(Node vertex) {
        CladePartitionObservation observation = new CladePartitionObservation(
                getBranchLengthOld(vertex),
                getBranchLengthOldOld(vertex),
                getBranchLengthOldYoung(vertex),
                getBranchLengthOldSmall(vertex),
                getBranchLengthOldBig(vertex),
                getHeightOld(vertex)
        );
        return observation;
    }

    @Override
    protected void increaseOccurrenceCount(Node vertex) {
        super.increaseOccurrenceCount(vertex);

        CladePartitionObservation observation = createObservation(vertex);
        this.observations.add(observation);
    }

    @Override
    protected void decreaseOccurrenceCount(Node vertex) {
        super.decreaseOccurrenceCount(vertex);

        CladePartitionObservation observation = createObservation(vertex);
        this.observations.removeIf(x -> x == observation);
    }

    public List<CladePartitionObservation> getObservations() {
        return observations;
    }

    /** -- Convenience getters - Convenience getters **/

    public DoubleStream getObservedBranchLengthsOld() {
        return this.getObservations().stream().mapToDouble(x -> x.branchLengthOld());
    }

    public DoubleStream getObservedLogBranchLengthsOld() {
        return this.getObservations().stream().mapToDouble(x -> Utils.logOrZero(x.branchLengthOld()));
    }

    /* -- DISTRIBUTION - DISTRIBUTION -- */

    protected Function<CladePartitionObservation, BranchLengthDistribution> distributionFunc;

    public Function<CladePartitionObservation, BranchLengthDistribution> getDistributionFunc() {
        return distributionFunc;
    }

    public void setDistributionFunc(Function<CladePartitionObservation, BranchLengthDistribution> distributionFunc) {
        this.distributionFunc = distributionFunc;
    }

    @Override
    public double getCCP(Node vertex) {
        double ccdCCP = super.getCCP();

        double minBranchLength = getBranchLengthOld(vertex);

        BranchLengthDistribution branchLengthDistribution = this.getBranchLengthDistribution(vertex);
        double branchProbability = branchLengthDistribution.density(minBranchLength);

        return ccdCCP * branchProbability;
    }

    @Override
    public double getLogCCP(Node vertex) {
        double logCcdCCP = super.getLogCCP();

        double minBranchLength = getBranchLengthOld(vertex);

        BranchLengthDistribution branchLengthDistribution = this.getBranchLengthDistribution(vertex);
        double logBranchProbability = branchLengthDistribution.logDensity(minBranchLength);

        return logCcdCCP + logBranchProbability;
    }

    /* -- Sampling - Sampling -- */

    protected BranchLengthDistribution getBranchLengthDistribution(Node vertex) {
        CladePartitionObservation observation = createObservation(vertex);
        return this.getDistributionFunc().apply(observation);
    }

    public double sampleMinBranchLength(Node vertex) {
        BranchLengthDistribution dist = this.getBranchLengthDistribution(vertex);
        return dist.sample();
    }

    /* -- MAP Tree - MAP Tree -- */
    public double getMAPBranchLength(Node vertex) {
        BranchLengthDistribution dist = this.getBranchLengthDistribution(vertex);
        try {
            return dist.mode();
        } catch (NoModeException e) {
            return dist.mean();
        }
    }
}
