package ccd.model;

import java.util.List;

public class BCCDMuSigmaLocalBetaOldOld extends BCCDLinear {
    public BCCDMuSigmaLocalBetaOldOld() {
        super(
                List.of(x -> Utils.logOrZero(x.branchLengthOldOld())),
                false
        );
    }

}
