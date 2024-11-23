package ccd.model;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

public class BetaDistribution extends BranchLengthDistribution {
    org.apache.commons.math3.distribution.BetaDistribution betaDistribution;

    public BetaDistribution(double alpha, double beta) {
        this.betaDistribution = new org.apache.commons.math3.distribution.BetaDistribution(alpha, beta);
    }

    @Override
    public double density(double value) {
        try {
            return this.betaDistribution.density(value);
        } catch (Exception e) {
            return 1.0;
        }
    }

    @Override
    public double mode() {
        double alpha = this.betaDistribution.getAlpha();
        double beta = this.betaDistribution.getBeta();

        if (alpha < 1 || beta < 1) {
            throw new AssertionError("Distribution has no mode.");
        }

        return (alpha - 1) / (alpha + beta - 2);
    }

    @Override
    public double sample() {
        return this.betaDistribution.sample();
    }

    public static BetaDistribution estimate(double[] observations) {
        if (observations.length == 1) return new BetaDistribution(1, 1);

        double mean = new Mean().evaluate(observations);
        double variance = new Variance().evaluate(observations);

        double alpha = mean * (mean * (1 - mean) / variance - 1);
        double beta = (1 - mean) * (mean * (1 - mean) / variance - 1);

        if (alpha < 1) alpha = 1.05;
        if (beta < 1) beta = 1.05;

        return new BetaDistribution(alpha, beta);
    }
}
