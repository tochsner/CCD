package ccd.model;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

public class BetaDistribution extends BranchLengthDistribution {
    org.apache.commons.math3.distribution.BetaDistribution betaDistribution;

    public BetaDistribution(double alpha, double beta) {
        this.betaDistribution = new org.apache.commons.math3.distribution.BetaDistribution(
                Math.max(0.3, alpha),
                Math.max(0.3, beta)
        );
    }

    @Override
    public double density(double value) {
        return this.betaDistribution.density(value);
    }

    @Override
    public double logDensity(double value) {
        return this.betaDistribution.logDensity(value);
    }

    @Override
    public double mode() throws NoModeException {
        double alpha = this.betaDistribution.getAlpha();
        double beta = this.betaDistribution.getBeta();

        if (alpha < 1 || beta < 1) {
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

        return mean * (mean * (1 - mean) / variance - 1);
    }

    public static double estimateBeta(double[] observations) {
        if (observations.length < 2)
            throw new IllegalArgumentException("Too few samples to estimate the distribution.");

        double mean = new Mean().evaluate(observations);
        double variance = new Variance().evaluate(observations);

        return (1 - mean) * (mean * (1 - mean) / variance - 1);
    }

    public static BetaDistribution estimate(double[] observations) {
        return new BetaDistribution(
                BetaDistribution.estimateAlpha(observations),
                BetaDistribution.estimateBeta(observations)
        );
    }
}
