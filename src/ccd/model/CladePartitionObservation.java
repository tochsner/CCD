package ccd.model;

public record CladePartitionObservation(
        double logBranchLengthOld,
        double logBranchLengthOldOld,
        double logBranchLengthOldYoung,
        double logBranchLengthOldSmall,
        double logBranchLengthOldBig
) { }
