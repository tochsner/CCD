package ccd.model;

import java.util.List;

public class BCCDMuSigmaLocalBetaOldSmallOldBig extends BCCDLinear {
    public BCCDMuSigmaLocalBetaOldSmallOldBig() {
        super(
                List.of(x -> x.getObservedLogBranchLengthsOldSmall(), x -> x.getObservedLogBranchLengthsOldBig()),
                List.of(x -> x.logBranchLengthOldSmall(), x -> x.logBranchLengthOldBig()),
                false
        );
    }
}
