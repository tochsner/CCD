package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Function;

public abstract class BCCDFeatureSelector {

    private final boolean logFeatures;

    public BCCDFeatureSelector() {
        this(false);
    }

    public BCCDFeatureSelector(boolean logFeatures) {
        this.logFeatures = logFeatures;
    }

    abstract double scoreFeature(double[] x, double[] y);


    public Function<CladePartitionObservation, Double>[] getFeatures(List<BCCDCladePartition> partitions) {
        Function<CladePartitionObservation, Double>[] features = new Function[partitions.size()];

        Function<CladePartitionObservation, Double>[] possibleFeatures;
        if (logFeatures)
            possibleFeatures = new Function[]{
                    (Function<CladePartitionObservation, Double>) v -> Utils.logOrZero(v.branchLengthOldOld()),
                    (Function<CladePartitionObservation, Double>) v -> Utils.logOrZero(v.branchLengthOldYoung()),
                    (Function<CladePartitionObservation, Double>) v -> Utils.logOrZero(v.branchLengthOldSmall()),
                    (Function<CladePartitionObservation, Double>) v -> Utils.logOrZero(v.branchLengthOldBig())
            };
        else
            possibleFeatures = new Function[]{
                    (Function<CladePartitionObservation, Double>) v -> v.branchLengthOldOld(),
                    (Function<CladePartitionObservation, Double>) v -> v.branchLengthOldYoung(),
                    (Function<CladePartitionObservation, Double>) v -> v.branchLengthOldSmall(),
                    (Function<CladePartitionObservation, Double>) v -> v.branchLengthOldBig()
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

            List<Double> scores = new ArrayList<>();
            for (double[] y : ys) {
                scores.add(this.scoreFeature(x, y));
            }

            if (Collections.max(scores) == Collections.min(scores)) {
                // all scores are the same
                idxWithoutEnoughData.add(i);
                continue;
            }

            double maxCov = Collections.max(scores);

            for (int j = 0; j < possibleFeatures.length; j++) {
                if (maxCov == scores.get(j)) {
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
