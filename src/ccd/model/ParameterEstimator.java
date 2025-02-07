package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.Collection;
import java.util.List;
import java.util.Map;

public abstract class ParameterEstimator<T extends AbstractCCD> {
    public abstract T buildCCD(int numLeaves, boolean storeBaseTrees);

    public abstract void estimateParameters(T ccd);

    public abstract int getNumberOfParameters(T ccd);

    public void estimateMAPBranches(Tree topology, List<BCCDCladePartition> partitions) {}

    public Map<Node, Double> sampleAllLengths(Map<Node, CTMCCCDCladePartition> partitions) {
        throw new UnsupportedOperationException();
    }

    public Map<Node, Double> getAllMAPLengths(Map<Node, CTMCCCDCladePartition> partitions) {
        throw new UnsupportedOperationException();
    }
}
