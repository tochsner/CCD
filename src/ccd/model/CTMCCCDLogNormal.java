package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.*;


public class CTMCCCDLogNormal extends ParameterEstimator<CTMCCCD> {

    @Override
    public CTMCCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new CTMCCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(CTMCCCD sbccd) {
        this.estimateTimeSinceLastSplitDistribution(sbccd);
        this.estimatePartitionRates(sbccd);
    }

    @Override
    public int getNumberOfParameters(CTMCCCD ccd) {
        return ccd.getNumberOfCladePartitions();
    }

    private void estimateTimeSinceLastSplitDistribution(CTMCCCD sbccd) {
        double[] timesSinceLastSplit = new double[sbccd.getNumberOfBaseTrees()];

        for (int i = 0; i < sbccd.getNumberOfBaseTrees(); i++) {
            Tree tree = sbccd.getBaseTrees().get(i);
            timesSinceLastSplit[i] = Math.log(Utils.getTimeSinceLastSplit(tree.getRoot()));
        }

        BranchLengthDistribution distribution = LogNormalDistribution.estimateMLE(
                timesSinceLastSplit
        );
        sbccd.setTimeSinceLastSplitDistribution(distribution);
    }


    private void estimatePartitionRates(CTMCCCD sbccd) {
        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Double> logMeans = new ArrayList<>();
        List<Double> logStds = new ArrayList<>();

        List<Clade> clades = sbccd.getClades().stream().toList();
        for (int i = 0; i < sbccd.getNumberOfClades(); i++) {
            CTMCCCDClade clade = (CTMCCCDClade) clades.get(i);

            if (clade.isRoot()) {
                continue;
            }

            if (clade.getNumberOfOccurrences() <= 5) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            double[] logObservations = clade.getPartitions().stream().flatMap(x -> ((CTMCCCDCladePartition) x).getObservations().stream()).mapToDouble(x -> Math.log(x.branchLengthTop())).toArray();
            LogNormalDistribution distribution = LogNormalDistribution.estimateMLE(logObservations);

            logMeans.add(distribution.logMean);
            logStds.add(distribution.logStd);

            for (CladePartition partition : clade.getPartitions()) {
                ((CTMCCCDCladePartition) partition).setTopBranchLengthDistribution(distribution);
            }
        }

        double logMean = logMeans.stream().mapToDouble(x -> x).average().orElseThrow();
        double logStd = logStds.stream().mapToDouble(x -> x).average().orElseThrow();

        for (int idx : idxWithoutEnoughData) {
            CTMCCCDClade clade = (CTMCCCDClade) clades.get(idx);
            LogNormalDistribution distribution = new LogNormalDistribution(logMean, logStd);

            for (CladePartition partition : clade.getPartitions()) {
                ((CTMCCCDCladePartition) partition).setTopBranchLengthDistribution(distribution);
            }
        }
    }

    @Override
    public Map<Node, Double> sampleAllLengths(Map<Node, CTMCCCDCladePartition> partitions) {
        Map<Node, Double> topBranchLengths = new HashMap<>();

        for (Node vertex : partitions.keySet()) {
            CTMCCCDCladePartition partition = partitions.get(vertex);

            double topBranchLength = partition.getParentClade().isRoot() ? 0.0 : partition.getTopBranchLengthDistribution().sample();
            topBranchLengths.put(vertex, topBranchLength);
        }

        return topBranchLengths;
    }

    @Override
    public Map<Node, Double> getAllMAPLengths(Map<Node, CTMCCCDCladePartition> partitions) {
        Map<Node, Double> topBranchLengths = new HashMap<>();

        for (Node vertex : partitions.keySet()) {
            CTMCCCDCladePartition partition = partitions.get(vertex);

            double topBranchLength;
            try {
                topBranchLength = partition.getParentClade().isRoot() ? 0.0 : partition.getTopBranchLengthDistribution().mode();
            } catch (NoModeException e) {
                topBranchLength = partition.getTopBranchLengthDistribution().mean();
            }
            topBranchLengths.put(vertex, topBranchLength);
        }

        return topBranchLengths;
    }

}
