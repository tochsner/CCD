package ccd.model;

import java.util.List;

public class BCCDMuSigmaBetaOldSmallOldBig extends BCCDLinear {
    public BCCDMuSigmaBetaOldSmallOldBig() {
        super(
                List.of(
                        x -> Utils.logOrZero(x.branchLengthOldSmall()),
                        x -> Utils.logOrZero(x.branchLengthOldBig())
                ),
                true
        );
    }
}
