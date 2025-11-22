package ccd.model;

import org.apache.commons.math3.analysis.function.Logit;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.DoubleStream;


public class SBCCDLogisticMultivariateGaussian extends ParameterEstimator<SBCCD> {

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
        Map<BitSet, Double> cladeMus = new HashMap<>();
        Map<BitSet, Double> cladeSigmas = new HashMap<>();

        for (Clade clade : sbccd.getClades()) {
            if (clade.isRoot() || clade.isLeaf()) continue;

            BitSet cladeBitSet = clade.getCladeInBits();

            double[] observedFractions = getObservedBranchFractions(clade);

            if (observedFractions.length <= 2) continue;

            LogitNormalDistribution dist = LogitNormalDistribution.estimateMLE(observedFractions);

            double cladeMu = dist.mean;
            double cladeSigma = dist.std;

            cladeMus.put(cladeBitSet, cladeMu);
            cladeSigmas.put(cladeBitSet, cladeSigma);
        }

        // set clade partition parameters

        double meanMu = cladeMus.values().stream().mapToDouble(x -> x).sum() / cladeMus.size();
        double meanSigma = cladeSigmas.values().stream().mapToDouble(x -> x).sum() / cladeSigmas.size();

        for (SBCCDCladePartition partition : sbccd.getAllPartitions()) {
            Clade firstClade = Utils.getFirstClade(partition);
            Clade secondClade = Utils.getSecondClade(partition);

            if (!firstClade.isLeaf()) {
                double mu = cladeMus.getOrDefault(firstClade.getCladeInBits(), meanMu);
                double sigma = cladeSigmas.getOrDefault(firstClade.getCladeInBits(), meanSigma);
                partition.setFirstBranchAlpha(mu);
                partition.setFirstBranchBeta(sigma);
            }

            if (!secondClade.isLeaf()) {
                double mu = cladeMus.getOrDefault(secondClade.getCladeInBits(), meanMu);
                double sigma = cladeSigmas.getOrDefault(secondClade.getCladeInBits(), meanSigma);
                partition.setSecondBranchAlpha(mu);
                partition.setSecondBranchBeta(sigma);
            }
        }
    }

    private double[] getObservedBranchFractions(Clade clade) {
        List<Double> observedFractions = new ArrayList<>();

        for (Clade parent : clade.getParentClades()) {
            for (CladePartition partition : parent.getPartitions()) {
                List<SBCCDCladePartitionObservation> observations = ((SBCCDCladePartition) partition).getObservations();
                for (SBCCDCladePartitionObservation observation : observations) {
                    if (Utils.getFirstClade(partition) == clade) {
                        observedFractions.add(observation.branchLengthLeft() / observation.subTreeHeight());
                    } else if (Utils.getSecondClade(partition) == clade) {
                        observedFractions.add(observation.branchLengthRight() / observation.subTreeHeight());
                    }
                }
            }
        }

        if (observedFractions.size() == 0) {
            throw new AssertionError("No parent clades found. This should not happen.");
        }

        Logit logit = new Logit();
        return observedFractions.stream().mapToDouble(x -> logit.value(x)).toArray();
    }
}
