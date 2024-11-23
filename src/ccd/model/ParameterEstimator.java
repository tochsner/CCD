package ccd.model;

import beast.base.evolution.tree.Tree;

import java.util.List;

public abstract class ParameterEstimator<T extends AbstractCCD> {
    public abstract T getCCD(int numLeaves, boolean storeBaseTrees);

    public abstract void estimateParameters(T ccd);

    public void estimateMAPBranches(Tree topology, List<BCCDCladePartition> partitions) {}
}
