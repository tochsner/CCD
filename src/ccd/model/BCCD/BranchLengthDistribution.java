package ccd.model.BCCD;

public abstract class BranchLengthDistribution {
    public abstract double sample();
    public abstract double density(double value);
    public abstract double mode();
}
