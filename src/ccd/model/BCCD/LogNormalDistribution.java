package ccd.model.BCCD;

import org.apache.commons.math3.distribution.NormalDistribution;
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
        double logNormalDensity = FastMath.exp(normalDensity);
        return logNormalDensity;
    }

    @Override
    public double mode() {
        return Math.exp(this.logMean);
    }

    @Override
    public double sample() {
        double sampledLogValue = new NormalDistribution(this.logMean, this.logStd).sample();
        return Math.exp(sampledLogValue);
    }
}
