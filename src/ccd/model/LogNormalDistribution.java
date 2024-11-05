package ccd.model;

import org.apache.commons.math3.distribution.NormalDistribution;

public class LogNormalDistribution extends NormalDistribution {
    public LogNormalDistribution(double logMean, double logStd) {
        super(logMean, logStd);
    }

    @Override
    public double density(double value) {
        double logValue = Math.log(value);
        return super.density(logValue);
    }

    @Override
    public double sample() {
        double sampledLogValue = super.sample();
        return Math.exp(sampledLogValue);
    }
}
