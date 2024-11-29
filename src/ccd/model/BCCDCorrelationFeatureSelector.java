package ccd.model;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class BCCDCorrelationFeatureSelector extends BCCDFeatureSelector {
    public BCCDCorrelationFeatureSelector() {
        super();
    }

    public BCCDCorrelationFeatureSelector(boolean logFeatures) {
        super(logFeatures);
    }

    @Override
    double scoreFeature(double[] x, double[] y) {
        PearsonsCorrelation correlation = new PearsonsCorrelation();
        return Math.abs(correlation.correlation(x, y));
    }
}
