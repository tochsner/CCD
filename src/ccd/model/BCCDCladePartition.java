package ccd.model;

import beast.base.evolution.tree.Node;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

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

    protected static double getLogBranchLength(Node vertex, int childIndex) {
        double vertexHeight = vertex.getHeight();
        double childHeight = vertex.getChild(childIndex).getHeight();

        if (vertexHeight - childHeight <= 0)
            throw new IllegalArgumentException();

        return Math.log(vertexHeight - childHeight);
    }

    protected static double getLogBranchLengthOld(Node vertex) {
        if (vertex.isLeaf()) return 0.0;
        int oldChild = getOldChildIndex(vertex);
        return getLogBranchLength(vertex, oldChild);
    }

    protected static double getLogBranchLengthYoung(Node vertex) {
        if (vertex.isLeaf()) return 0.0;
        int youngChild = 1 - getOldChildIndex(vertex);
        return getLogBranchLength(vertex, youngChild);
    }

    protected static double getLogBranchLengthOldOld(Node vertex) {
        if (vertex.getLeft().isLeaf() || vertex.getRight().isLeaf()) return 0.0;
        Node oldChild = vertex.getChild(getOldChildIndex(vertex));
        return getLogBranchLengthOld(oldChild);
    }

    protected static double getLogBranchLengthOldYoung(Node vertex) {
        if (vertex.getLeft().isLeaf() || vertex.getRight().isLeaf()) return 0.0;
        Node oldChild = vertex.getChild(getOldChildIndex(vertex));
        return getLogBranchLengthYoung(oldChild);
    }

    protected static double getLogBranchLengthOldSmall(Node vertex) {
        if (vertex.getLeft().isLeaf() || vertex.getRight().isLeaf()) return 0.0;
        Node oldChild = vertex.getChild(getOldChildIndex(vertex));
        if (oldChild.isLeaf()) return 0.0;
        int smallChild = getSmallChildIndex(oldChild);
        return getLogBranchLength(oldChild, smallChild);
    }

    protected static double getLogBranchLengthOldBig(Node vertex) {
        if (vertex.getLeft().isLeaf() || vertex.getRight().isLeaf()) return 0.0;
        Node oldChild = vertex.getChild(getOldChildIndex(vertex));
        if (oldChild.isLeaf()) return 0.0;
        int bigChild = 1 - getSmallChildIndex(oldChild);
        return getLogBranchLength(oldChild, bigChild);
    }

    private static CladePartitionObservation createObservation(Node vertex) {
        CladePartitionObservation observation = new CladePartitionObservation(
                getLogBranchLengthOld(vertex),
                getLogBranchLengthOldOld(vertex),
                getLogBranchLengthOldYoung(vertex),
                getLogBranchLengthOldSmall(vertex),
                getLogBranchLengthOldBig(vertex)
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

    public DoubleStream getObservedLogBranchLengthsOld() {
        return this.getObservations().stream().mapToDouble(x -> x.logBranchLengthOld());
    }

    public DoubleStream getObservedLogBranchLengthsOldOld() {
        return this.getObservations().stream().mapToDouble(x -> x.logBranchLengthOldOld());
    }

    public DoubleStream getObservedLogBranchLengthsOldYoung() {
        return this.getObservations().stream().mapToDouble(x -> x.logBranchLengthOldYoung());
    }

    public DoubleStream getObservedLogBranchLengthsOldSmall() {
        return this.getObservations().stream().mapToDouble(x -> x.logBranchLengthOldSmall());
    }

    public DoubleStream getObservedLogBranchLengthsOldBig() {
        return this.getObservations().stream().mapToDouble(x -> x.logBranchLengthOldBig());
    }

    /* -- DISTRIBUTION - DISTRIBUTION -- */

    protected Function<CladePartitionObservation, Double> logMeanFunc;
    protected Function<CladePartitionObservation, Double> logVarianceFunc;

    public Function<CladePartitionObservation, Double> getLogMeanFunc() {
        return logMeanFunc;
    }

    public void setLogMeanFunc(Function<CladePartitionObservation, Double> logMeanFunc) {
        this.logMeanFunc = logMeanFunc;
    }

    public Function<CladePartitionObservation, Double> getLogVarianceFunc() {
        return logVarianceFunc;
    }

    public void setLogVarianceFunc(Function<CladePartitionObservation, Double> logVarianceFunc) {
        this.logVarianceFunc = logVarianceFunc;
    }

    @Override
    public double getCCP(Node vertex) {
        double ccdCCP = super.getCCP();

        double minBranchLength = getLogBranchLengthOld(vertex);

        AbstractRealDistribution branchLengthDistribution = this.getBranchLengthDistribution(vertex);
        double branchProbability = branchLengthDistribution.density(minBranchLength);

        return ccdCCP * branchProbability;
    }

    /* -- Sampling - Sampling -- */

    protected AbstractRealDistribution getBranchLengthDistribution(Node vertex) {
        CladePartitionObservation observation = createObservation(vertex);

        double logMean = this.getLogMeanFunc().apply(observation);
        double logVariance = this.getLogVarianceFunc().apply(observation);

        return new NormalDistribution(logMean, Math.sqrt(logVariance));
    }

    public double sampleMinBranchLength(Node vertex) {
        AbstractRealDistribution dist = this.getBranchLengthDistribution(vertex);
        return Math.exp(dist.sample());
    }
}
