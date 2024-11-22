package ccd.model;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.IllinoisSolver;
import org.apache.commons.math3.exception.NoBracketingException;

public class InverseTrigamma {
    public static double value(double y) {
        UnivariateFunction trigammaFunction = v -> Trigamma.value(v) - y;

        double initialGuess = 1 / y;

        IllinoisSolver solver = new IllinoisSolver();
        try {
            return solver.solve(100, trigammaFunction, initialGuess - 100, initialGuess + 100);
        } catch (NoBracketingException e) {
            throw new ArithmeticException("No solution found in the interval for value " + y);
        }
    }
}
