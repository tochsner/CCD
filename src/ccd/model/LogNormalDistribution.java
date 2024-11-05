package ccd.model;

import org.apache.commons.math3.distribution.NormalDistribution;

public class LogNormalDistribution extends BranchLengthDistribution {
    NormalDistribution normalDistribution;

    public LogNormalDistribution(double logMean, double logStd) {
        this.normalDistribution = new NormalDistribution(logMean, logStd);
    }

    @Override
    public double density(double value) {
        double logValue = Math.log(value);
        return this.normalDistribution.density(logValue);
    }

    @Override
    public double mode() {
        return Math.exp(this.normalDistribution.getMean());
    }

    @Override
    public double sample() {
        double sampledLogValue = this.normalDistribution.sample();
        return Math.exp(sampledLogValue);
    }
}
