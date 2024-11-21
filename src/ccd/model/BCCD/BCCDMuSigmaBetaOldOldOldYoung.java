package ccd.model.BCCD;

import java.util.List;

public class BCCDMuSigmaBetaOldOldOldYoung extends BCCDLinear {
    public BCCDMuSigmaBetaOldOldOldYoung() {
        super(
                List.of(
                        x -> Utils.logOrZero(x.branchLengthOldOld()),
                        x -> Utils.logOrZero(x.branchLengthOldYoung())
                ),
                true
        );
    }
}
