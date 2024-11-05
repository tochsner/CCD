package ccd.model;

import java.util.ArrayList;
import java.util.List;

/**
 * Not finished yet.
 */
public class BCCDGammaMuSigmaBetaOldOld extends BCCDGammaLinear {
    public BCCDGammaMuSigmaBetaOldOld() {
        super(List.of(x -> x.branchLengthOldOld()), true);
    }
}
