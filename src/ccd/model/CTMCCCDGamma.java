package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class CTMCCCDGamma extends ParameterEstimator<CTMCCCD> {

    @Override
    public CTMCCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        return new CTMCCCD(numLeaves, true, this);
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
            timesSinceLastSplit[i] = Utils.getTimeSinceLastSplit(tree.getRoot());
        }

        BranchLengthDistribution distribution = GammaDistribution.estimateMLE(
                timesSinceLastSplit, 1e-6
        );
        sbccd.setTimeSinceLastSplitDistribution(distribution);
    }

    private void estimatePartitionRates(CTMCCCD sbccd) {
        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Double> shapes = new ArrayList<>();
        List<Double> scales = new ArrayList<>();

        List<Clade> clades = sbccd.getClades().stream().toList();
        for (int i = 0; i < sbccd.getNumberOfClades(); i++) {
            CTMCCCDClade clade = (CTMCCCDClade) clades.get(i);

            if (clade.isRoot()) {
                continue;
            }
            if (clade.isLeaf()) {
                continue;
            }

            if (clade.getNumberOfOccurrences() <= 5) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            double[] observations = clade.getPartitions().stream().flatMap(x -> ((CTMCCCDCladePartition) x).getObservations().stream()).mapToDouble(x -> x.branchLengthTop()).toArray();
            GammaDistribution distribution = GammaDistribution.estimateMLE(observations, 1e-5);

            shapes.add(distribution.gammaDistribution.getShape());
            scales.add(distribution.gammaDistribution.getScale());

            for (CladePartition partition : clade.getPartitions()) {
                ((CTMCCCDCladePartition) partition).setTopBranchLengthDistribution(distribution);
            }
        }

        double shape = shapes.stream().mapToDouble(x -> x).average().orElseThrow();
        double scale = scales.stream().mapToDouble(x -> x).average().orElseThrow();

        for (int idx : idxWithoutEnoughData) {
            CTMCCCDClade clade = (CTMCCCDClade) clades.get(idx);
            GammaDistribution distribution = new GammaDistribution(shape, scale);

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

            double topBranchLength = partition.getParentClade().isRoot() ? 0.0 : partition.getTopBranchLengthDistribution().mean();
            topBranchLengths.put(vertex, topBranchLength);
        }

        return topBranchLengths;
    }

}
