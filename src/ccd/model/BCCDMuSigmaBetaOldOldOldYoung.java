package ccd.model;

import java.util.List;

public class BCCDMuSigmaBetaOldOldOldYoung extends BCCDLinear {
    public BCCDMuSigmaBetaOldOldOldYoung() {
        super(
                List.of(x -> x.getObservedLogBranchLengthsOldOld(), x -> x.getObservedLogBranchLengthsOldYoung()),
                List.of(x -> x.logBranchLengthOldOld(), x -> x.logBranchLengthOldYoung())
        );
    }
}
