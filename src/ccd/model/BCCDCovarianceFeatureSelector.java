package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Function;

public class BCCDCovarianceFeatureSelector extends BCCDFeatureSelector {
    @Override
    public Function<CladePartitionObservation, Double>[] getFeatures(List<BCCDCladePartition> partitions) {
        Function<CladePartitionObservation, Double>[] features = new Function[partitions.size()];

        Function<CladePartitionObservation, Double>[] possibleFeatures = new Function[]{
                (Function<CladePartitionObservation, Double>) v -> Utils.logOrZero(v.branchLengthOldOld()),
                (Function<CladePartitionObservation, Double>) v -> Utils.logOrZero(v.branchLengthOldYoung()),
                (Function<CladePartitionObservation, Double>) v -> Utils.logOrZero(v.branchLengthOldSmall()),
                (Function<CladePartitionObservation, Double>) v -> Utils.logOrZero(v.branchLengthOldBig())
        };

        List<Integer> idxWithoutEnoughData = new LinkedList<>();
        List<Integer> numPicks = new ArrayList<>();

        for (Function feature : possibleFeatures)
            numPicks.add(0);

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            if (partition.getNumberOfOccurrences() == 1) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            double[] x = partition.getObservations().stream().mapToDouble(v -> v.branchLengthOld()).toArray();

            List<double[]> ys = new ArrayList<>();
            for (int j = 0; j < possibleFeatures.length; j++) {
                Function<CladePartitionObservation, Double> feature = possibleFeatures[j];
                ys.add(partition.getObservations().stream().mapToDouble(v -> feature.apply(v)).toArray());
            }

            Covariance covariance = new Covariance();
            List<Double> absCovariances = new ArrayList<>();
            for (double[] y : ys) {
                absCovariances.add(Math.abs(covariance.covariance(x, y)));
            }

            if (Collections.max(absCovariances) == Collections.min(absCovariances)) {
                // all absolute covariances are the same
                idxWithoutEnoughData.add(i);
                continue;
            }

            double maxCov = Collections.max(absCovariances);

            for (int j = 0; j < possibleFeatures.length; j++) {
                if (maxCov == absCovariances.get(j)) {
                    numPicks.set(j, numPicks.get(j) + 1);
                    features[i] = possibleFeatures[j];
                    break;
                }
            }
        }

        double maxPicks = Collections.max(numPicks);
        Function<CladePartitionObservation, Double> mostCommonFunc = null;

        for (int j = 0; j < possibleFeatures.length; j++) {
            if (maxPicks == numPicks.get(j)) {
                mostCommonFunc = possibleFeatures[j];
                break;
            }
        }

        for (int idx : idxWithoutEnoughData) {
            features[idx] = mostCommonFunc;
        }

        return features;
    }
}
