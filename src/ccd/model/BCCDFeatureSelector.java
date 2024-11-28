package ccd.model;

import java.util.List;
import java.util.function.Function;

public abstract class BCCDFeatureSelector {
    public abstract Function<CladePartitionObservation, Double>[] getFeatures(List<BCCDCladePartition> partitions);
}
