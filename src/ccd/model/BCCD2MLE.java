package ccd.model;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;

import java.io.Serializable;
import java.util.List;
import java.util.Map;

/**
 * This class implements the log MLE for the BCCD2 model (up to some constant).
 * <p>
 * Format of the parameter vector: [mu_s ... sigma_s ... beta]
 */
public class BCCD2MLE {

    List<BCCD2CladePartition> partitions;

    public BCCD2MLE(List<BCCD2CladePartition> partitions) {
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

    public ObjectiveFunctionGradient logMLEGradient() {
        return new ObjectiveFunctionGradient(parameters -> {
            double beta = parameters[parameters.length - 1];

            double[] gradient = new double[parameters.length];
            for (int i = 0; i < this.partitions.size(); i++) {
                BCCD2CladePartition partition = this.partitions.get(i);

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
}
