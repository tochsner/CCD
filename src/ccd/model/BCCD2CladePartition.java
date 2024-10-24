package ccd.model;

import beast.base.evolution.tree.Node;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ConstantRealDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;

import java.util.HashMap;
import java.util.Map;

public class BCCD2CladePartition extends BCCDCladePartition {
    public BCCD2CladePartition(BCCDClade parentClade, BCCDClade[] childClades) {
        super(parentClade, childClades);
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
        Node firstChild = vertex.getChild(0);
        Node secondChild = vertex.getChild(1);

        double minBranchLengthDown;
        if (firstChild.isLeaf() || secondChild.isLeaf()) {
            minBranchLengthDown = 0;
        } else {
            int minBranchChild = getMinBranchChildIdx(vertex);
            minBranchLengthDown = getMinBranchLength(vertex.getChild(minBranchChild));
        }

        return getMu() + getBeta() * Math.log(minBranchLengthDown);
    }

    @Override
    protected double getLogVariance(Node vertex, double logSampleMean) {
        return this.getSigma();
    }
}
