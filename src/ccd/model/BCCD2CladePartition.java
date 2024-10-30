package ccd.model;

import beast.base.evolution.tree.Node;

import java.util.*;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

public class BCCD2CladePartition extends BCCDCladePartition {
    private List<CladePartitionObservation> observations = new LinkedList<>();

    public BCCD2CladePartition(BCCD2Clade parentClade, BCCD2Clade[] childClades) {
        super(parentClade, childClades);
    }

    /* -- BOOK KEEPING - BOOK KEEPING -- */

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

    /* -- DISTRIBUTION PARAMETERS - DISTRIBUTION PARAMETERS -- */

    protected double mu;
    protected double sigma;
    protected double beta;

    public double getMu() {
        return mu;
    }

    public void setMu(double mu) {
        this.mu = mu;
    }

    public double getSigma() {
        return sigma;
    }

    public void setSigma(double sigma) {
        this.sigma = sigma;
    }

    public double getBeta() {
        return beta;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }

    /* -- Sampling - Sampling -- */

    @Override
    protected double getLogMean(Node vertex) {
        return getMu() + getBeta() * getMinLogBranchLengthDown(vertex);
    }

    @Override
    protected double getLogVariance(Node vertex, double logSampleMean) {
        return this.getSigma();
    }
}
