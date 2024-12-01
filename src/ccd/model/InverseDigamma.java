package ccd.model;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.IllinoisSolver;
import org.apache.commons.math3.exception.NoBracketingException;

public class InverseDigamma {
    public static double value(double y) {
        UnivariateFunction digammaFunction = v -> Digamma.value(v) - y;

        // see https://arxiv.org/pdf/1705.06547
        double initialGuess = 1 / Math.log(1 + Math.exp(-y));

        IllinoisSolver solver = new IllinoisSolver(1e-5);
        try {
            return solver.solve(
                    500,
                    digammaFunction,
                    initialGuess - 0.5,
                    initialGuess + 0.5,
                    initialGuess
            );
        } catch (NoBracketingException e) {
            return initialGuess;
            // throw new ArithmeticException("No solution found in the interval for value " + y);
        }
    }
}
