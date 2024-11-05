package ccd.model;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.util.List;

public class BCCDMLE extends BCCDParameterEstimator {
    @Override
    public void estimateParameters(List<BCCDCladePartition> partitions) {
        NonLinearConjugateGradientOptimizer optimizer = new NonLinearConjugateGradientOptimizer(
                NonLinearConjugateGradientOptimizer.Formula.FLETCHER_REEVES,
                new SimpleValueChecker(1e-4, 0)
        );

        double[] solution = optimizer.optimize(
                GoalType.MAXIMIZE,
                new InitialGuess(this.getInitialGuess(partitions)),
                this.logMLE(partitions),
                this.logMLEGradient(partitions),
                new MaxEval(50000)
        ).getPoint();

        this.setPartitionParameters(partitions, solution);
    }

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
                    double b = observation.logBranchLengthOld();
                    double bDown = observation.logBranchLengthOldOld();

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
                    double b = observation.logBranchLengthOld();
                    double bDown = observation.logBranchLengthOldOld();

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

    public double[] getInitialGuess(List<BCCDCladePartition> partitions) {
        double[] initialParameters = new double[2*partitions.size() + 1];

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double[] branchLengths = partition.getObservedLogBranchLengthsOld().toArray();
            initialParameters[i] = new Mean().evaluate(branchLengths);
            initialParameters[partitions.size() + i] = new Variance().evaluate(branchLengths);
        }

        initialParameters[initialParameters.length - 1] = 0.0;

        return initialParameters;
    }

    protected static void setPartitionParameters(List<BCCDCladePartition> partitions, double[] solution) {
        double beta = solution[solution.length - 1];

        for (int i = 0; i < partitions.size(); i++) {
            double mu = solution[i];
            double sigma = solution[partitions.size() + i];

            BCCDCladePartition partition = partitions.get(i);
            partition.setDistributionFunc(x -> new NormalDistribution(mu + beta * x.logBranchLengthOldOld(), Math.sqrt(sigma)));
        }
    }
}
