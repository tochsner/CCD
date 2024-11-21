package ccd.model.BCCD;

import java.util.List;

public class BCCDMuSigmaBetaOldOld extends BCCDLinear {
    public BCCDMuSigmaBetaOldOld() {
        super(
                List.of(x -> Utils.logOrZero(x.branchLengthOldOld())),
                true
        );
    }

}
