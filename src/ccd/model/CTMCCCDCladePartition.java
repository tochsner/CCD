package ccd.model;

import beast.base.evolution.tree.Node;

import java.util.LinkedList;
import java.util.List;

public class CTMCCCDCladePartition extends CladePartition {
    private List<CTMCCCDCladePartitionObservation> observations = new LinkedList<>();

    public CTMCCCDCladePartition(CTMCCCDClade parentClade, CTMCCCDClade[] childClades) {
        super(parentClade, childClades);
    }

    /* -- BOOK KEEPING - BOOK KEEPING -- */

    protected static double getBranchLength(Node parent, Node child) {
        if (child.isRoot()) return 0.0;

        double vertexHeight = parent.getHeight();
        double childHeight = child.getHeight();

        if (vertexHeight - childHeight < 0)
            throw new IllegalArgumentException("Negative branch length found.");

        return vertexHeight - childHeight;
    }

    public static CTMCCCDCladePartitionObservation createObservation(Node vertex) {
        return new CTMCCCDCladePartitionObservation(
                getBranchLength(vertex.getParent(), vertex)
        );
    }

    @Override
    protected void increaseOccurrenceCount(Node vertex) {
        super.increaseOccurrenceCount(vertex);
        CTMCCCDCladePartitionObservation observation = createObservation(vertex);
        this.observations.add(observation);
    }

    @Override
    protected void decreaseOccurrenceCount(Node vertex) {
        super.decreaseOccurrenceCount(vertex);
        CTMCCCDCladePartitionObservation observation = createObservation(vertex);
        this.observations.removeIf(x -> x == observation);
    }

    public List<CTMCCCDCladePartitionObservation> getObservations() {
        return observations;
    }

    /* -- DISTRIBUTION - DISTRIBUTION -- */

    BranchLengthDistribution topBranchLengthDistribution;

    public BranchLengthDistribution getTopBranchLengthDistribution() {
        return topBranchLengthDistribution;
    }

    public void setTopBranchLengthDistribution(BranchLengthDistribution topBranchLengthDistribution) {
        this.topBranchLengthDistribution = topBranchLengthDistribution;
    }

    public double getCCP(Node vertex) {
        double ccdCCP = super.getCCP(vertex);

        if (vertex.isRoot()) {
            return ccdCCP;
        }

        CTMCCCDCladePartitionObservation observation = createObservation(vertex);
        BranchLengthDistribution branchLengthDistribution = this.getTopBranchLengthDistribution();
        double branchLengthProbability = branchLengthDistribution.density(observation.branchLengthTop());

        return ccdCCP + branchLengthProbability;
    }

    public double getLogCCP(Node vertex) {
        double logCcdCCP = super.getLogCCP(vertex);

        if (vertex.isRoot()) {
            return logCcdCCP;
        }

        CTMCCCDCladePartitionObservation observation = createObservation(vertex);
        BranchLengthDistribution branchLengthDistribution = this.getTopBranchLengthDistribution();
        double logBranchLengthProbability = branchLengthDistribution.logDensity(observation.branchLengthTop());

        return logCcdCCP + logBranchLengthProbability;
    }
}
