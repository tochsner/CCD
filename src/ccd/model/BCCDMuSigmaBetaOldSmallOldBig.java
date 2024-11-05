package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class BCCDMuSigmaBetaOldSmallOldBig extends BCCDLinear {
    public BCCDMuSigmaBetaOldSmallOldBig() {
        super(
                List.of(x -> x.getObservedLogBranchLengthsOldSmall(), x -> x.getObservedLogBranchLengthsOldBig()),
                List.of(x -> x.logBranchLengthOldSmall(), x -> x.logBranchLengthOldBig()),
                true
        );
    }
}
