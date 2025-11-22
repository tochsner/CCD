package ccd.model;

import org.apache.commons.math3.analysis.function.Sigmoid;
import org.apache.commons.math3.analysis.function.Logit;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

public class LogitNormalDistribution extends BranchLengthDistribution {
    double mean;
    double std;

    Logit logit = new Logit();
    Sigmoid logistic = new Sigmoid();

    public LogitNormalDistribution(double mean, double std) {
        this.mean = mean;
        this.std = std;
    }

    @Override
    public double density(double value) {
        return new org.apache.commons.math3.distribution.NormalDistribution(this.mean, this.std).density(value);
    }

    @Override
    public double logDensity(double value) {
        return new org.apache.commons.math3.distribution.NormalDistribution(this.mean, this.std).logDensity(value);
    }

    @Override
    public double mode() {
        return logistic.value(this.mean);
    }

    @Override
    public double mean() {
        return logistic.value(this.mean);
    }

    @Override
    public double sample() {
        return logistic.value(this.mean + this.std * this.getRandom().nextGaussian());
    }

    public static LogitNormalDistribution estimateMLE(double[] observations) {
        double mean = new Mean().evaluate(observations);
        double variance = new Variance().evaluate(observations);
        return new LogitNormalDistribution(mean, Math.sqrt(variance));
    }
}
