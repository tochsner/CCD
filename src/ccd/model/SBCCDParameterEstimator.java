package ccd.model;

import beast.base.evolution.tree.Tree;

import java.util.List;

public abstract class SBCCDParameterEstimator {
    public abstract void estimateParameters(SBCCD bccd);

    public void estimateMAPBranches(Tree topology, List<SBCCDCladePartition> partitions) {}
}
