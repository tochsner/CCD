package ccd.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.DoubleStream;


public class SBCCDCladeBeta extends ParameterEstimator<SBCCD> {

    double K;

    public SBCCDCladeBeta() {
        this(10.0);
    }

    public SBCCDCladeBeta(double K) {
        this.K = K;
    }

    @Override
    public SBCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new SBCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(SBCCD sbccd) {
        this.estimateHeightDistribution(sbccd);
        this.estimatePartitionParameters(sbccd);
    }

    @Override
    public int getNumberOfParameters(SBCCD ccd) {
        int heightParameters = 2;

        int fractionParameters = 0;
        for (Clade clade : ccd.getClades()) {
            if (!clade.isLeaf()) fractionParameters++;
        }

        return heightParameters + fractionParameters;
    }

    private void estimateHeightDistribution(SBCCD sbccd) {
        Clade rootClade = sbccd.getRootClade();

        DoubleStream observedHeights = DoubleStream.empty();

        for (SBCCDCladePartition partition : sbccd.getAllPartitions()) {
            if (partition.getParentClade() == rootClade) {
                observedHeights = DoubleStream.concat(
                        observedHeights,
                        partition.getObservations().stream().mapToDouble(x -> Utils.logOrZero(x.subTreeHeight()))
                );
            }
        }

        double[] observedLogHeightsArray = observedHeights.toArray();
        LogNormalDistribution heightDistribution = LogNormalDistribution.estimateMLE(observedLogHeightsArray);
        sbccd.setHeightDistribution(heightDistribution);
    }

    private void estimatePartitionParameters(SBCCD sbccd) {
        Map<BitSet, Double> cladeAlphas = new HashMap<>();
        Map<BitSet, Double> cladeBetas = new HashMap<>();

        for (Clade clade : sbccd.getClades()) {
            if (clade.isRoot()) continue;

            BitSet cladeBitSet = clade.getCladeInBits();

            double[] observedFractions = getObservedBranchFractions(clade);

            if (observedFractions.length <= 2) continue;

            double cladeAlpha = BetaDistribution.estimateAlpha(observedFractions);
            double cladeBeta = BetaDistribution.estimateBeta(observedFractions);

            cladeAlphas.put(cladeBitSet, cladeAlpha);
            cladeBetas.put(cladeBitSet, cladeBeta);
        }

        assert cladeAlphas.size() == sbccd.getNumberOfClades();

        // set clade partition parameters

        double meanAlpha = cladeAlphas.values().stream().mapToDouble(x -> x).sum() / cladeAlphas.size();
        double meanBeta = cladeBetas.values().stream().mapToDouble(x -> x).sum() / cladeAlphas.size();

        for (SBCCDCladePartition partition : sbccd.getAllPartitions()) {
            Clade firstClade = partition.getChildClades()[0];
            Clade secondClade = partition.getChildClades()[1];

            if (!firstClade.isLeaf()) {
                double alpha = cladeAlphas.getOrDefault(firstClade.getCladeInBits(), meanAlpha);
                double beta = cladeBetas.getOrDefault(firstClade.getCladeInBits(), meanBeta);
                partition.setFirstBranchAlpha(alpha);
                partition.setFirstBranchBeta(beta);
            }

            if (!secondClade.isLeaf()) {
                double alpha = cladeAlphas.getOrDefault(secondClade.getCladeInBits(), meanAlpha);
                double beta = cladeBetas.getOrDefault(secondClade.getCladeInBits(), meanBeta);
                partition.setSecondBranchAlpha(alpha);
                partition.setSecondBranchBeta(beta);
            }
        }
    }

    private double[] getObservedBranchFractions(Clade clade) {
        List<Double> observedFractions = new ArrayList<>();

        for (Clade parent : clade.getParentClades()) {
            for (CladePartition partition : parent.getPartitions()) {
                List<SBCCDCladePartitionObservation> observations = ((SBCCDCladePartition) partition).getObservations();
                for (SBCCDCladePartitionObservation observation : observations) {
                    if (partition.getChildClades()[0] == clade) {
                        observedFractions.add(observation.branchLengthLeft() / observation.subTreeHeight());
                    } else if (partition.getChildClades()[1] == clade) {
                        observedFractions.add(observation.branchLengthRight() / observation.subTreeHeight());
                    }
                }
            }
        }

        if (observedFractions.size() == 0) {
            throw new AssertionError("No parent clades found. This should not happen.");
        }

        return observedFractions.stream().mapToDouble(x -> x).toArray();
    }
}
