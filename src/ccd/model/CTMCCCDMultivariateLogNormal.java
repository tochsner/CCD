package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.*;


public class CTMCCCDMultivariateLogNormal extends ParameterEstimator<CTMCCCD> {

    private List<Clade> clades;

    Map<CladePartition, Integer> partitionIdx = new HashMap<>();
    double[] logMeans;
    double[][] logStd;

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
            timesSinceLastSplit[i] = Math.log(Utils.getTimeSinceLastSplit(tree.getRoot()));
        }

        BranchLengthDistribution distribution = LogNormalDistribution.estimateMLE(
                timesSinceLastSplit
        );
        sbccd.setTimeSinceLastSplitDistribution(distribution);
    }

    private void estimatePartitionRates(CTMCCCD sbccd) {
        this.clades = sbccd.getClades().stream().map(x -> x).toList();

        RealMatrix X = new BlockRealMatrix(sbccd.getNumberOfBaseTrees(), sbccd.getNumberOfClades());
        RealMatrix occurrences = new BlockRealMatrix(sbccd.getNumberOfBaseTrees(), sbccd.getNumberOfClades());

        for (int i = 0; i < sbccd.getNumberOfBaseTrees(); i++) {
            Tree tree = sbccd.getBaseTrees().get(i);

            Map<Node, CTMCCCDClade> collectedClades = new HashMap<>();
            sbccd.collectClades(tree.getRoot(), collectedClades);

            for (Node vertex : collectedClades.keySet()) {
                CTMCCCDClade clade = collectedClades.get(vertex);

                CTMCCCDCladePartitionObservation observation = CTMCCCDCladePartition.createObservation(vertex);
                double topBranchLength = observation.branchLengthTop();

                X.setEntry(i, this.clades.indexOf(clade), Math.log(topBranchLength));
                occurrences.setEntry(i, this.clades.indexOf(clade), 1.0);
            }
        }

        this.logMeans = new double[sbccd.getNumberOfClades()];

        for (int i = 0; i < sbccd.getNumberOfClades(); i++) {
            int count = 0;
            double sum = 0.0;

            for (int j = 0; j < sbccd.getNumberOfBaseTrees(); j++) {
                if (occurrences.getEntry(j, i) == 0.0) continue;
                count++;
                sum += X.getEntry(j, i);
            }

            this.logMeans[i] = sum / count;
        }

        this.logStd = new double[sbccd.getNumberOfClades()][sbccd.getNumberOfClades()];

        for (int i = 0; i < sbccd.getNumberOfClades(); i++) {
            for (int j = 0; j < sbccd.getNumberOfClades(); j++) {
                int count = 0;
                double sum = 0.0;

                for (int t = 0; t < sbccd.getNumberOfBaseTrees(); t++) {
                    if (occurrences.getEntry(t, i) == 0.0 || occurrences.getEntry(t, j) == 0.0) continue;
                    count++;
                    sum += (X.getEntry(t, i) - this.logMeans[i]) * (X.getEntry(t, j) - this.logMeans[j]);
                }

                this.logStd[i][j] = sum / count;
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

            int i = this.clades.indexOf(partition.getParentClade());
            double topBranchLength = partition.getParentClade().isRoot() ? 0.0 : Math.exp(this.logMeans[i]);
            topBranchLengths.put(vertex, topBranchLength);
        }

        return topBranchLengths;
    }

}
