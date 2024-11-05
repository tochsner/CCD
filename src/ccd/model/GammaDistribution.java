package ccd.model;

import org.apache.commons.math3.distribution.NormalDistribution;

public class GammaDistribution extends BranchLengthDistribution {
    org.apache.commons.math3.distribution.GammaDistribution gammaDistribution;

    public GammaDistribution(double shape, double scale) {
        this.gammaDistribution = new org.apache.commons.math3.distribution.GammaDistribution(shape, scale);
    }

    @Override
    public double density(double value) {
        return this.gammaDistribution.density(value);
    }

    @Override
    public double mode() {
        double shape = this.gammaDistribution.getShape();
        double scale = this.gammaDistribution.getScale();

        if (shape < 1) {
            return 0.0;
        } else {
            return (shape - 1) * scale;
        }
    }

    @Override
    public double sample() {
        return this.gammaDistribution.sample();
    }
}
