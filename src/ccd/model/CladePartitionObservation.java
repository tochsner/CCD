package ccd.model;

public record CladePartitionObservation(
        double logBranchLengthOld,
        double logBranchLengthYoung,
        double logBranchLengthOldOld,
        double logBranchLengthOldYoung,
        double logBranchLengthOldSmall,
        double logBranchLengthOldBig
) { }
