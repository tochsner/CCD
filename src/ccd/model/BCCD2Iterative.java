package ccd.model;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;

import java.util.List;
import java.util.Map;

/**
 * This class implements the log MLE for the BCCD2 model (up to some constant).
 * <p>
 * Format of the parameter vector: [mu_s ... sigma_s ... beta]
 */
public class BCCD2Iterative extends BCCD2MLE {

    public BCCD2Iterative(List<BCCD2CladePartition> partitions) {
        super(partitions);
    }

    public void updateMusSigmas(double[] parameters) {
        double beta = parameters[parameters.length - 1];

        int a = 0;

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);
            int numObservations = partition.getObservations().size();

            double mu = 0;
            for (CladePartitionObservation observation : partition.getObservations()) {
                double b = observation.logMinBranchLength();
                double bDown = observation.logMinBranchLengthDown();

                mu += b - beta*bDown;

                a++;
            }
            mu /= numObservations;

            double sigma = 0;
            for (CladePartitionObservation observation : partition.getObservations()) {
                double b = observation.logMinBranchLength();
                double bDown = observation.logMinBranchLengthDown();

                sigma += Math.pow(b - beta*bDown - mu, 2);
            }
            sigma /= numObservations;

            parameters[i] = mu;
            parameters[i + this.partitions.size()] = sigma;
        }
    }

    public void updateBeta(double[] parameters) {
        double nominator = 0;
        double denominator = 0;

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);

            double mu = parameters[i];
            double sigma = parameters[partitions.size() + i];

            if (partition.getObservations().size() == 1) continue;

            for (CladePartitionObservation observation : partition.getObservations()) {
                double b = observation.logMinBranchLength();
                double bDown = observation.logMinBranchLengthDown();

                nominator += (b - mu) * bDown / sigma;
                denominator += Math.pow(bDown, 2) / sigma;
            }
        }

        parameters[parameters.length - 1] = nominator / denominator;
    }

    public double[] estimateParameters() {
        SimpleValueChecker convergenceChecker = new SimpleValueChecker(1e-6, 0);

        double[] solution = new double[2 * this.partitions.size() + 1];

        solution[solution.length - 1] = 0.1;

        int iteration = 0;
        double previousMLE;
        double currentMLE = this.logMLE().getObjectiveFunction().value(solution);

        while (true) {
            this.updateMusSigmas(solution);
            this.updateBeta(solution);

            previousMLE = currentMLE;
            currentMLE = this.logMLE().getObjectiveFunction().value(solution);

            boolean hasConverged = convergenceChecker.converged(
                    iteration++,
                    new PointValuePair(solution, previousMLE),
                    new PointValuePair(solution, currentMLE)
            );

            if (hasConverged)
                break;
        }

        return solution;
    }
}
