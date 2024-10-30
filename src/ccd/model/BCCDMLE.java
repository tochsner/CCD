package ccd.model;

import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;

import java.util.List;

public class BCCDMLE extends BCCDParameterEstimator {

    public ObjectiveFunction logMLE(List<BCCDCladePartition> partitions) {
        return new ObjectiveFunction(parameters -> {
            double beta = parameters[parameters.length - 1];

            double logMLE = 0.0;
            for (int i = 0; i < partitions.size(); i++) {
                BCCDCladePartition partition = partitions.get(i);

                double mu = parameters[i];
                double sigma = parameters[partitions.size() + i];

                if (sigma == 0.0) continue;

                for (CladePartitionObservation observation : partition.getObservations()) {
                    double b = observation.logMinBranchLength();
                    double bDown = observation.logMinBranchLengthDown();

                    logMLE += 0.5 * (-Math.log(2*sigma*Math.PI) - Math.pow(b - mu - beta * bDown, 2) / sigma);
                }
            }

            return logMLE;
        });
    }

    public ObjectiveFunctionGradient logMLEGradient(List<BCCDCladePartition> partitions) {
        return new ObjectiveFunctionGradient(parameters -> {
            double beta = parameters[parameters.length - 1];

            double[] gradient = new double[parameters.length];
            for (int i = 0; i < partitions.size(); i++) {
                BCCDCladePartition partition = partitions.get(i);

                double mu = parameters[i];
                double sigma = parameters[partitions.size() + i];

                if (sigma == 0.0) continue;

                for (CladePartitionObservation observation : partition.getObservations()) {
                    double b = observation.logMinBranchLength();
                    double bDown = observation.logMinBranchLengthDown();

                    // mu
                    gradient[i] += (b - mu - beta * bDown) / sigma;

                    // sigma
                    gradient[partitions.size() + i] += 0.5 * (
                            (-1 / sigma) + Math.pow(b - mu - beta * bDown, 2) / Math.pow(sigma, 2)
                    );

                    // beta
                    gradient[parameters.length - 1] += (b - mu - beta * bDown) * bDown / sigma;
                }
            }
            return gradient;
        });
    }

    @Override
    public void estimateParameters(List<BCCDCladePartition> partitions) {

    }
}
