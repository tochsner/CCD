package ccd.model;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.FastMath;

public class ExponentialDistribution extends BranchLengthDistribution {
    org.apache.commons.math3.distribution.ExponentialDistribution distribution;

    public ExponentialDistribution(double mean) {
        this.distribution = new org.apache.commons.math3.distribution.ExponentialDistribution(mean);
    }

    @Override
    public double density(double value) {
        return this.distribution.density(value);
    }

    @Override
    public double logDensity(double value) {
        return this.distribution.logDensity(value);
    }

    @Override
    public double mode() throws NoModeException {
        throw new NoModeException();
    }

    @Override
    public double mean() {
        return this.distribution.getMean();
    }

    @Override
    public double sample() {
       return this.distribution.sample();
    }

    public static ExponentialDistribution estimateMLE(double[] observations) {
        double mean = new Mean().evaluate(observations);
        return new ExponentialDistribution(mean);
    }
}
