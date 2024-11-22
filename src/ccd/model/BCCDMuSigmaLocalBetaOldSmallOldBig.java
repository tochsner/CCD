package ccd.model;

import java.util.List;

public class BCCDMuSigmaLocalBetaOldSmallOldBig extends BCCDLinear {
    public BCCDMuSigmaLocalBetaOldSmallOldBig() {
        super(
                List.of(
                        x -> Utils.logOrZero(x.branchLengthOldSmall()),
                        x -> Utils.logOrZero(x.branchLengthOldBig())
                ),
                false
        );
    }
}
