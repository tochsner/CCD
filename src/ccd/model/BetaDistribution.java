package ccd.model;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

public class BetaDistribution extends BranchLengthDistribution {
    org.apache.commons.math3.distribution.BetaDistribution betaDistribution;

    public BetaDistribution(double alpha, double beta) {
        this.betaDistribution = new org.apache.commons.math3.distribution.BetaDistribution(
                Math.max(1e-5, alpha),
                Math.max(1e-5, beta)
        );
    }

    @Override
    public double density(double value) {
        return this.betaDistribution.density(
                Math.max(1e-5, Math.min(1 - 1e-5, value))
        );
    }

    @Override
    public double logDensity(double value) {
        return this.betaDistribution.logDensity(
                Math.max(1e-5, Math.min(1 - 1e-5, value))
        );
    }

    @Override
    public double mode() throws NoModeException {
        double alpha = this.betaDistribution.getAlpha();
        double beta = this.betaDistribution.getBeta();

        if (alpha <= 1 || beta <= 1) {
            throw new NoModeException();
        }

        return (alpha - 1) / (alpha + beta - 2);
    }

    @Override
    public double mean() {
        return this.betaDistribution.getNumericalMean();
    }

    @Override
    public double sample() {
        return this.betaDistribution.sample();
    }

    public static double estimateAlpha(double[] observations) {
        if (observations.length < 2)
            throw new IllegalArgumentException("Too few samples to estimate the distribution.");

        double mean = new Mean().evaluate(observations);
        double variance = new Variance().evaluate(observations);

        if (variance >= mean * (1 - mean))
            throw new IllegalArgumentException("Beta parameters could not be estimated.");

        return mean * (mean * (1 - mean) / variance - 1);
    }

    public static double estimateBeta(double[] observations) {
        if (observations.length < 2)
            throw new IllegalArgumentException("Too few samples to estimate the distribution.");

        double mean = new Mean().evaluate(observations);
        double variance = new Variance().evaluate(observations);

        if (variance >= mean * (1 - mean))
            throw new IllegalArgumentException("Beta parameters could not be estimated.");

        return (1 - mean) * (mean * (1 - mean) / variance - 1);
    }

    public static BetaDistribution estimate(double[] observations) {
        return new BetaDistribution(
                BetaDistribution.estimateAlpha(observations),
                BetaDistribution.estimateBeta(observations)
        );
    }

    public static double[] estimateMAP(double[] observations, double relTolerance) {
        PowellOptimizer optimizer = new PowellOptimizer(
                relTolerance, relTolerance, new SimpleValueChecker(relTolerance, relTolerance)
        );

        double initialAlpha;
        double initialBeta;

        if (observations.length < 2) {
            initialAlpha = 2.0;
            initialBeta = 2.0;
        } else {
            initialAlpha = BetaDistribution.estimateAlpha(observations);
            initialBeta = BetaDistribution.estimateBeta(observations);
        }

        InitialGuess initialGuess = new InitialGuess(new double[]{
                FastMath.log(initialAlpha), FastMath.log(initialBeta)
        });

        PointValuePair result = optimizer.optimize(
                new MaxEval(100000),
                initialGuess,
                GoalType.MAXIMIZE,
                initialGuess,
                new ObjectiveFunction(new BetaWithPrior(observations))
        );

        return new double[]{
                FastMath.exp(result.getPoint()[0]),
                FastMath.exp(result.getPoint()[1])
        };
    }

    private static class BetaWithPrior implements MultivariateFunction {
        double[] observations;

        public BetaWithPrior(double[] observations) {
            this.observations = observations;
        }

        public double value(double[] variables) {
            double alpha = Math.exp(variables[0]);
            double beta = Math.exp(variables[1]);

            if (alpha == 0.0 || beta == 0.0) return Float.NEGATIVE_INFINITY;
            if (1000 < alpha || 1000 < beta) return Float.NEGATIVE_INFINITY;

            double z = Gamma.logGamma(alpha) + Gamma.logGamma(beta) - Gamma.logGamma(alpha + beta);

            double logP = 0;

            for (double observation : this.observations) {
                logP += (alpha - 1) * FastMath.log(observation) + (beta - 1) * FastMath.log(1 - observation) - z;
            }

            // mu = alpha / (alpha + beta)
            // we have a uniform prior on mu, so we don't have to add anything

            double logVar = FastMath.log(alpha * beta / FastMath.pow(alpha + beta, 2) / (alpha + beta + 1));

            // we have a lognormal(-2, 1) prior on phi
            double priorLogP = -FastMath.pow(logVar + 2, 2) / 2 - logVar - FastMath.log(FastMath.sqrt(2 * FastMath.PI));
            logP += priorLogP;

            return logP;
        }
    }
}
