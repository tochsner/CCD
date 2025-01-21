package ccd.model;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.IllinoisSolver;
import org.apache.commons.math3.exception.NoBracketingException;

import java.util.*;
import java.util.stream.DoubleStream;


public class SBCCDDirichlet extends ParameterEstimator<SBCCD> {

    Double K = null;

    public SBCCDDirichlet() {
    }

    public SBCCDDirichlet(double K) {
        this.K = K;
    }

    @Override
    public SBCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new SBCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(SBCCD sbccd) {
        Double K = this.K;
        if (K == null) {
            K = estimateK(sbccd);
        }

        this.estimateHeightDistribution(sbccd);
        this.estimatePartitionParameters(sbccd, K);
    }

    private static double estimateK(SBCCD sbccd) {
        double sumK = 0.0;
        double normalization = 0.0;

        for (SBCCDCladePartition partition : sbccd.getAllPartitions()) {
            if (partition.getNumberOfOccurrences() < 5) {
                continue;
            }

            if (!Utils.getFirstClade(partition).isLeaf()) {
                double[] firstBranchRatios = partition.getObservations().stream().mapToDouble(x -> x.branchLengthLeft() / x.subTreeHeight()).toArray();

                double alpha = BetaDistribution.estimateAlpha(firstBranchRatios);
                double beta = BetaDistribution.estimateBeta(firstBranchRatios);

                double K = alpha + beta;

                if (!Double.isNaN(K)) {
                    sumK += K * partition.getNumberOfOccurrences();
                    normalization += partition.getNumberOfOccurrences();
                }
            }

            if (!Utils.getSecondClade(partition).isLeaf()) {
                double[] secondBranchRatios = partition.getObservations().stream().mapToDouble(x -> x.branchLengthRight() / x.subTreeHeight()).toArray();

                double alpha = BetaDistribution.estimateAlpha(secondBranchRatios);
                double beta = BetaDistribution.estimateBeta(secondBranchRatios);

                double K = alpha + beta;

                if (!Double.isNaN(K)) {
                    sumK += K * partition.getNumberOfOccurrences();
                    normalization += partition.getNumberOfOccurrences();
                }
            }
        }

        return sumK / normalization;
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

    private static void estimateHeightDistribution(SBCCD sbccd) {
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

    private static void estimatePartitionParameters(SBCCD sbccd, double K) {
        Map<BitSet, Double> cladeAlphas = new HashMap<>();

        for (Clade outerClade : sbccd.getClades()) {
            for (Clade clade : sbccd.getClades()) {
                BitSet cladeBitSet = clade.getCladeInBits();

                if (cladeAlphas.containsKey(cladeBitSet)) continue;

                boolean hasUnknownParentClade = false;

                for (Clade parentClade : clade.getParentClades()) {
                    if (!cladeAlphas.containsKey(parentClade.getCladeInBits())) {
                        hasUnknownParentClade = true;
                        break;
                    }
                }

                if (hasUnknownParentClade) continue;

                double cladeAlpha = estimateCladeAlpha(clade, cladeAlphas, K);
                cladeAlphas.put(cladeBitSet, cladeAlpha);
            }
        }

        assert cladeAlphas.size() == sbccd.getNumberOfClades();

        // set clade partition alphas

        for (SBCCDCladePartition partition : sbccd.getAllPartitions()) {
            BitSet parentClade = partition.getParentClade().getCladeInBits();

            Clade firstClade = Utils.getFirstClade(partition);
            Clade secondClade = Utils.getSecondClade(partition);

            double parentAlpha = cladeAlphas.get(parentClade);

            if (!firstClade.isLeaf()) {
                double firstAlpha = cladeAlphas.get(firstClade.getCladeInBits());
                partition.setFirstBranchAlpha(parentAlpha - firstAlpha);
                partition.setFirstBranchBeta(firstAlpha);
            }

            if (!secondClade.isLeaf()) {
                double secondAlpha = cladeAlphas.get(secondClade.getCladeInBits());
                partition.setSecondBranchAlpha(parentAlpha - secondAlpha);
                partition.setSecondBranchBeta(secondAlpha);
            }
        }
    }

    private static double estimateCladeAlpha(Clade clade, Map<BitSet, Double> existingCladeAlphas, double K) {
        if (clade.isRoot()) {
            return K;
        }

        if (clade.isLeaf()) {
            return 0.0;
        }

        if (clade.getNumberOfOccurrences() < 2) {
            return existingCladeAlphas.values().stream().mapToDouble(x -> x).average().orElseThrow();
        }

        List<Double> parentAlphas = new ArrayList<>();
        List<Double> observedFractions = new ArrayList<>();

        for (Clade parent : clade.getParentClades()) {
            double parentAlpha = existingCladeAlphas.get(parent.getCladeInBits());

            for (CladePartition partition : parent.getPartitions()) {
                List<SBCCDCladePartitionObservation> observations = ((SBCCDCladePartition) partition).getObservations();
                for (SBCCDCladePartitionObservation observation : observations) {
                    if (Utils.getFirstClade(partition) == clade) {
                        observedFractions.add(observation.branchLengthLeft() / observation.subTreeHeight());
                        parentAlphas.add(parentAlpha);
                    } else if (Utils.getSecondClade(partition) == clade) {
                        observedFractions.add(observation.branchLengthRight() / observation.subTreeHeight());
                        parentAlphas.add(parentAlpha);
                    }
                }
            }
        }

        if (observedFractions.size() == 0) {
            throw new AssertionError("No parent clades found. This should not happen.");
        }

        int n = observedFractions.size();
        double logF = observedFractions.stream().mapToDouble(x -> Math.log(x)).sum();
        double logOneMinusF = observedFractions.stream().mapToDouble(x -> Math.log(1 - x)).sum();

        UnivariateFunction alphaFunction = alpha -> n * Digamma.value(alpha) - parentAlphas.stream().mapToDouble(x -> Digamma.value(x - alpha)).sum() + logF - logOneMinusF;

        double maxAlpha = parentAlphas.stream().mapToDouble(x -> x).min().orElseThrow();
        double initialGuess = maxAlpha / 2;
        IllinoisSolver solver = new IllinoisSolver(1e-5);
        try {
            double alpha = solver.solve(
                    500,
                    alphaFunction,
                    1e-4,
                    maxAlpha - 1e-4,
                    initialGuess
            );
            return alpha;
        } catch (NoBracketingException e) {
            System.err.println("No solution found for Dirichlet parameter.");
            return maxAlpha - 1e-4;
        }
    }
}
