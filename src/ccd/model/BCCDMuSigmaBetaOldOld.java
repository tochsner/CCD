package ccd.model;

import java.util.List;

public class BCCDMuSigmaBetaOldOld extends BCCDLinear {
    public BCCDMuSigmaBetaOldOld() {
        super(
                List.of(x -> x.getObservedLogBranchLengthsOldOld()),
                List.of(x -> x.logBranchLengthOldOld())
        );
    }

}
