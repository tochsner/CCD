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

    private static CTMCCCDCladePartitionObservation createObservation(Node vertex) {
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

    protected double rate;

    public void setRate(double rate) {
        this.rate = rate;
    }

    public double getRate() {
        return this.rate;
    }

    public double getCCP(Node vertex) {
        throw new UnsupportedOperationException();
    }

    public double getLogCCP(Node vertex) {
        if (vertex.isRoot()) {
            return Math.log(1.0 * this.getNumberOfOccurrences() / this.getParentClade().getNumberOfOccurrences());
        }

        double sumRates = 0.0;
        for (CladePartition partition : this.getParentClade().getPartitions()) {
            sumRates += ((CTMCCCDCladePartition) partition).getRate();
        }
        double logTransitionProbability = Math.log(this.getRate() / sumRates);

        CTMCCCDCladePartitionObservation observation = createObservation(vertex);
        double branchLengthTop = observation.branchLengthTop();
        double logDwellTimeProbability = Math.log(sumRates) - sumRates * branchLengthTop;

        return logTransitionProbability + logDwellTimeProbability;
    }
}
