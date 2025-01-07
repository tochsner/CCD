package ccd.model;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;


public class CCD1Estimator extends ParameterEstimator<CCD1> {
    @Override
    public CCD1 buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new CCD1(numLeaves, storeBaseTrees);
    }

    @Override
    public void estimateParameters(CCD1 ccd) {
        // nothing to do here
    }

    @Override
    public int getNumberOfParameters(CCD1 ccd) {
        return (int) ccd.getNumberOfParameters();
    }
}
