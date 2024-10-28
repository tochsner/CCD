package ccd.model;

import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;

import java.util.List;
import java.util.Map;

/**
 * This class implements the log MLE for the BCCD2 model (up to some constant).
 * <p>
 * Format of the parameter vector: [mu_s ... sigma_s ... beta]
 */
public class BCCD2Iterative {

    List<BCCD2CladePartition> partitions;

    public BCCD2Iterative(List<BCCD2CladePartition> partitions) {
        this.partitions = partitions;
    }

    public double[] getMusSigmas(double beta) {
        double[] initialParameters = new double[2*this.partitions.size() + 1];

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);

            initialParameters[i] = partition.getLogMeanApproximation();
            initialParameters[partitions.size() + i] = 2*partition.getLogVarianceApproximation();
        }

        return initialParameters;
    }


    public double[] getInitialGuess(double beta) {
        double[] initialParameters = new double[2*this.partitions.size() + 1];

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);

            initialParameters[i] = 2*partition.getLogMeanApproximation();
            initialParameters[partitions.size() + i] = 2*partition.getLogVarianceApproximation();
        }

        initialParameters[initialParameters.length - 1] = beta;

        return initialParameters;
    }

    public double getInitialBeta(double[] parameters) {
        double nominator = 0;
        double denominator = 0;

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);

            double mu = parameters[i];
            double sigma = parameters[partitions.size() + i];

            if (partition.getObservations().size() == 1) continue;

            for (Map.Entry<CladePartitionObservation, Integer> observationEntry : partition.getObservations().entrySet()) {
                int n = observationEntry.getValue();
                double b = observationEntry.getKey().logMinBranchLength();
                double bDown = observationEntry.getKey().logMinBranchLengthDown();

                nominator += (b - mu) * bDown / sigma;
                denominator +=Math.pow(bDown, 2) / sigma;
            }
        }

        return nominator / denominator;
    }
}
