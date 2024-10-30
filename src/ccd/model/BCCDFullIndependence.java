package ccd.model;

import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.util.List;

public class BCCDFullIndependence extends BCCDParameterEstimator {
    @Override
    public void estimateParameters(List<BCCDCladePartition> partitions) {
        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double mu = partition.getObservedMinLogBranchLength().average().getAsDouble();
            double sigma = new Variance().evaluate(partition.getObservedMinLogBranchLength().toArray());

            partition.setLogMeanFunc(x -> mu);
            partition.setLogVarianceFunc(x -> sigma);
        }
    }
}
