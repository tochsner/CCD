package ccd.model;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.MersenneTwister;

public abstract class BranchLengthDistribution {
    static RandomGenerator random;

    RandomGenerator getRandom() {
        if (BranchLengthDistribution.random == null) {
            BranchLengthDistribution.random = new MersenneTwister();
        }
        return BranchLengthDistribution.random;
    }

    public abstract double sample();
    public abstract double density(double value);
    public abstract double mode();
}
