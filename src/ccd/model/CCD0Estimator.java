package ccd.model;


public class CCD0Estimator extends ParameterEstimator<CCD0> {
    @Override
    public CCD0 buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new CCD0(numLeaves, storeBaseTrees);
    }

    @Override
    public void estimateParameters(CCD0 ccd) {
        // nothing to do here
    }

    @Override
    public int getNumberOfParameters(CCD0 ccd) {
        return (int) ccd.getNumberOfParameters();
    }
}
