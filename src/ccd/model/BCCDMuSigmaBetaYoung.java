package ccd.model;

import java.util.List;

public class BCCDMuSigmaBetaYoung extends BCCDLinear {
    public BCCDMuSigmaBetaYoung() {
        super(
                List.of(x -> x.getObservedLogBranchLengthsYoung()),
                List.of(x -> x.logBranchLengthYoung())
        );
    }
}
