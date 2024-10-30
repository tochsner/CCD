package ccd.model;

import java.util.List;

public abstract class BCCDParameterEstimator {
    public abstract void estimateParameters(List<BCCDCladePartition> partitions);
}
