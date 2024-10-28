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

    public ObjectiveFunction logMLE() {
        return new ObjectiveFunction(parameters -> {
            double beta = parameters[parameters.length - 1];

            double logMLE = 0.0;
            for (int i = 0; i < this.partitions.size(); i++) {
                BCCD2CladePartition partition = this.partitions.get(i);

                double mu = parameters[i];
                double sigma = parameters[partitions.size() + i];

                if (partition.getObservations().size() == 1) continue;

                for (Map.Entry<CladePartitionObservation, Integer> observationEntry : partition.getObservations().entrySet()) {
                    int n = observationEntry.getValue();
                    double b = observationEntry.getKey().logMinBranchLength();
                    double bDown = observationEntry.getKey().logMinBranchLengthDown();

                    logMLE += n * 0.5 * (-Math.log(sigma) - Math.pow(b - mu - beta * bDown, 2) / sigma);
                }
            }

            return logMLE;
        });
    }

    public ObjectiveFunctionGradient logMLEGradient() {
        return new ObjectiveFunctionGradient(parameters -> {
            double beta = parameters[parameters.length - 1];

            double[] gradient = new double[parameters.length];
            for (int i = 0; i < this.partitions.size(); i++) {
                BCCD2CladePartition partition = this.partitions.get(i);

                double mu = parameters[i];
                double sigma = parameters[partitions.size() + i];

                if (partition.getObservations().size() == 1) continue;

                for (Map.Entry<CladePartitionObservation, Integer> observationEntry : partition.getObservations().entrySet()) {
                    int n = observationEntry.getValue();
                    double b = observationEntry.getKey().logMinBranchLength();
                    double bDown = observationEntry.getKey().logMinBranchLengthDown();

                    // mu
                    gradient[i] += n * (b - mu - beta * bDown) / sigma;

                    // sigma
                    gradient[partitions.size() + i] += n * 0.5 * (
                            (-1 / sigma) + Math.pow(b - mu - beta * bDown, 2) / Math.pow(sigma, 2)
                    );

                    // beta
//                    gradient[parameters.length - 1] += n * (b - mu - beta * bDown) * bDown / sigma;
                }
            }
            return gradient;
        });
    }

    public double[] getInitialGuess() {
        double[] initialParameters = new double[2*this.partitions.size() + 1];

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);

            initialParameters[i] = 2*partition.getLogMeanApproximation();
            initialParameters[partitions.size() + i] = 2*partition.getLogVarianceApproximation();
        }

        initialParameters[initialParameters.length - 1] = 0.0;

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