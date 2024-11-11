package ccd.model;

public record CladePartitionObservation(
        double branchLengthOld,
        double branchLengthOldOld,
        double branchLengthOldYoung,
        double branchLengthOldSmall,
        double branchLengthOldBig,
        double heightOld
) { }
