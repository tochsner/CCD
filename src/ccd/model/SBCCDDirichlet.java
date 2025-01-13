package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.IllinoisSolver;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.sql.SQLOutput;
import java.util.*;
import java.util.stream.DoubleStream;


public class SBCCDDirichlet extends ParameterEstimator<SBCCD> {

    double K;

    public SBCCDDirichlet() {
        this(10.0);
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
        this.estimateHeightDistribution(sbccd);
        this.estimatePartitionParameters(sbccd);
    }

    @Override
    public int getNumberOfParameters(SBCCD ccd) {
        int heightParameters = 2;

        int fractionParameters = 0;
        for (Clade clade : ccd.getClades()) {
            if (!clade.isLeaf()) fractionParameters ++;
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

                double cladeAlpha = this.estimateCladeAlpha(clade, cladeAlphas);
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

    private double estimateCladeAlpha(Clade clade, Map<BitSet, Double> existingCladeAlphas) {
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

        UnivariateFunction alphaFunction = alpha -> n*Digamma.value(alpha) - parentAlphas.stream().mapToDouble(x -> Digamma.value(x - alpha)).sum() + logF - logOneMinusF;

        double maxAlpha = parentAlphas.stream().mapToDouble(x -> x).min().orElseThrow();
        double initialGuess = maxAlpha / 2;
        IllinoisSolver solver = new IllinoisSolver(1e-5);
        try {
            double alpha = solver.solve(
                    500,
                    alphaFunction,
                    0.01,
                    maxAlpha - 0.01,
                    initialGuess
            );
            return alpha;
        } catch (NoBracketingException e) {
            return maxAlpha - 0.1;
//            throw new ArithmeticException("No solution found in the interval for value.");
        }
    }
}
