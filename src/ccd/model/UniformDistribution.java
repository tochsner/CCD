package ccd.model;

public class UniformDistribution extends BranchLengthDistribution {
    @Override
    public double density(double value) {
       return 1.0;
    }

    @Override
    public double mode() {
        throw new UnsupportedOperationException();
    }

    @Override
    public double sample() {
        throw new UnsupportedOperationException();
    }
}
