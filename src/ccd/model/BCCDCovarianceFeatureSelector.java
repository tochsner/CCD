package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;

public class BCCDCovarianceFeatureSelector extends BCCDFeatureSelector {

    public BCCDCovarianceFeatureSelector() {
        super();
    }

    public BCCDCovarianceFeatureSelector(boolean logFeatures) {
        super(logFeatures);
    }

    @Override
    double scoreFeature(double[] x, double[] y) {
        Covariance covariance = new Covariance();
        return Math.abs(covariance.covariance(x, y));
    }
}
