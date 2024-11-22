package ccd.model;

import beast.base.evolution.tree.Tree;

import java.util.List;

public abstract class BCCDParameterEstimator {
    public abstract void estimateParameters(BCCD bccd);

    public void estimateMAPBranches(Tree topology, List<BCCDCladePartition> partitions) {}
}
