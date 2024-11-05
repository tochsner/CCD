package ccd.model;

import java.util.List;

public class BCCDMuSigmaLocalBetaOldOldOldYoung extends BCCDLinear {
    public BCCDMuSigmaLocalBetaOldOldOldYoung() {
        super(
                List.of(x -> x.getObservedLogBranchLengthsOldOld(), x -> x.getObservedLogBranchLengthsOldYoung()),
                List.of(x -> x.logBranchLengthOldOld(), x -> x.logBranchLengthOldYoung()),
                false
        );
    }
}
