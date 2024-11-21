package ccd.model.BCCD;

import java.util.List;

public class BCCDMuSigmaLocalBetaOldOldOldYoung extends BCCDLinear {
    public BCCDMuSigmaLocalBetaOldOldOldYoung() {
        super(
                List.of(
                        x -> Utils.logOrZero(x.branchLengthOldOld()),
                        x -> Utils.logOrZero(x.branchLengthOldYoung())
                ),
                false
        );
    }
}
