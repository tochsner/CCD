package ccd.model;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import java.util.Arrays;

public class GammaDistribution extends BranchLengthDistribution {
    org.apache.commons.math3.distribution.GammaDistribution gammaDistribution;

    public GammaDistribution(double shape, double scale) {
        this.gammaDistribution = new org.apache.commons.math3.distribution.GammaDistribution(
                this.getRandom(),
                shape,
                scale,
                1.0E-9
        );
    }

    @Override
    public double density(double value) {
        return this.gammaDistribution.density(value);
    }

    @Override
    public double logDensity(double value) {
        return this.gammaDistribution.logDensity(value);
    }

    @Override
    public double mode() throws NoModeException {
        double shape = this.gammaDistribution.getShape();
        double scale = this.gammaDistribution.getScale();

        if (shape < 1) {
            return 0;
        } else {
            return (shape - 1) * scale;
        }
    }

    @Override
    public double mean() {
        return this.gammaDistribution.getNumericalMean();
    }

    @Override
    public double sample() {
        return this.gammaDistribution.sample();
    }

    public static GammaDistribution estimateMLE(double[] observations, double relTolerance) {
        double shape = GammaDistribution.estimateShapeMLE(observations, relTolerance);
        double scale = GammaDistribution.estimateScaleMLE(observations, shape);
        return new GammaDistribution(shape, scale);
    }

    public static double estimateShapeMLE(double[] observations, double relTolerance) {
        double mean = new Mean().evaluate(observations);
        double logMean = Arrays.stream(observations).map(x -> Math.log(x)).average().orElseThrow();

        double oldShape = 0;
        double newShape = 0.5 / (Math.log(mean) - logMean);

        int MAX_ITERATIONS = 1000;
        int i = 0;
        while (relTolerance < Math.abs(oldShape - newShape) / oldShape) {
            oldShape = newShape;
            newShape = 1 / (
                    (1 / oldShape) + (logMean - Math.log(mean) + Math.log(oldShape) - Digamma.value(oldShape)) / (Math.pow(oldShape, 2) * (1 / oldShape - Trigamma.value(oldShape)))
            );

            i++;
            if (MAX_ITERATIONS <= i) {
                throw new ArithmeticException("Gamma estimation did not converge.");
            }
        }

        return newShape;
    }

    public static double estimateScaleMLE(double[] observations, double shapeMLE) {
        double mean = new Mean().evaluate(observations);
        return mean / shapeMLE;
    }
}
