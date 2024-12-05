package ccd.model;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.FastMath;

public class LogNormalDistribution extends BranchLengthDistribution {
    double logMean;
    double logStd;

    public LogNormalDistribution(double logMean, double logStd) {
        this.logMean = logMean;
        this.logStd = logStd;
    }

    @Override
    public double density(double value) {
        double logValue = FastMath.log(value);

        double x0 = logValue - this.logMean;
        double x1 = x0 / this.logStd;
        double logDensity = -0.5 * x1 * x1 - FastMath.log(this.logStd) + 0.5 * FastMath.log(6.283185307179586);
        double normalDensity = FastMath.exp(logDensity);
        double logNormalDensity = normalDensity / value;

        return logNormalDensity;
    }

    @Override
    public double logDensity(double value) {
        double logValue = FastMath.log(value);

        double x0 = logValue - this.logMean;
        double x1 = x0 / this.logStd;
        double logDensity = -0.5 * x1 * x1 - FastMath.log(this.logStd) + 0.5 * FastMath.log(6.283185307179586);
        double logNormalDensity = logDensity - logValue;

        return logNormalDensity;
    }

    @Override
    public double mode() {
        return Math.exp(this.logMean - this.logStd*this.logStd);
    }

    @Override
    public double mean() {
        return Math.exp(this.logMean);
    }

    @Override
    public double sample() {
        double sampledLogValue = this.logMean + this.logStd * this.getRandom().nextGaussian();
        return Math.exp(sampledLogValue);
    }

    public static LogNormalDistribution estimateMLE(double[] observations) {
        double logMean = new Mean().evaluate(observations);
        double logVariance = new Variance().evaluate(observations);
        return new LogNormalDistribution(logMean, Math.sqrt(logVariance));
    }
}
