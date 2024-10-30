package ccd.model;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleValueChecker;

import java.util.List;

public class BCCDIterativeMLE extends BCCDMLE {

    public void updateMusSigmas(List<BCCDCladePartition> partitions, double[] parameters) {
        double beta = parameters[parameters.length - 1];

        int a = 0;

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);
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
            parameters[i + partitions.size()] = sigma;
        }
    }

    public void updateBeta(List<BCCDCladePartition> partitions, double[] parameters) {
        double nominator = 0;
        double denominator = 0;

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

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

    @Override
    public void estimateParameters(List<BCCDCladePartition> partitions) {
        SimpleValueChecker convergenceChecker = new SimpleValueChecker(1e-6, 0);

        double[] solution = new double[2 * partitions.size() + 1];

        solution[solution.length - 1] = 0.1;

        int iteration = 0;
        double previousMLE;
        double currentMLE = this.logMLE(partitions).getObjectiveFunction().value(solution);

        while (true) {
            this.updateMusSigmas(partitions, solution);
            this.updateBeta(partitions, solution);

            previousMLE = currentMLE;
            currentMLE = this.logMLE(partitions).getObjectiveFunction().value(solution);

            boolean hasConverged = convergenceChecker.converged(
                    iteration++,
                    new PointValuePair(solution, previousMLE),
                    new PointValuePair(solution, currentMLE)
            );

            if (hasConverged)
                break;
        }

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);
            partition.setLogMean(solution[i]);
            partition.setLogVariance(solution[partitions.size() + i]);
        }
    }
}
