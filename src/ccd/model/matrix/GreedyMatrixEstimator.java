package ccd.model.matrix;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeUtils;
import ccd.model.HeightSettingStrategy;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.Graph;
import org.jgrapht.alg.interfaces.SpanningTreeAlgorithm;
import org.jgrapht.alg.spanning.PrimMinimumSpanningTree;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.SimpleWeightedGraph;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class GreedyMatrixEstimator extends MatrixEstimator {
    @Override
    public TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd) {
        Tree mapTopology = ccd.getMAPTree(HeightSettingStrategy.CommonAncestorHeights);
        TreeMatrixConfiguration configuration = this.getGreedyConfiguration(mapTopology, ccd.getBaseTrees());
        return configuration;
    }

    TreeMatrixConfiguration getGreedyConfiguration(Tree tree, List<Tree> observedTrees) {
        int n = tree.getLeafNodeCount();

        double[][] pairwiseScoresMatrix = new double[n][n];
        double[][][] pairwiseScoresPerTree = new double[n][n][observedTrees.size()];
        Node[][][] mrcas = new Node[n][n][observedTrees.size()];

        for (int t = 0; t < observedTrees.size(); t++) {
            Tree observedTree = observedTrees.get(t);
            SimpleGraph<Integer, DefaultEdge> treeGraph = CubeUtils.createUnweightedGraphForTree(observedTree);
            int[][] treeMrcas = CubeUtils.getMrcas(treeGraph, n, observedTree.getRoot().getNr());

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == j) continue;

                    Node mrca = observedTree.getNode(treeMrcas[i][j]);
                    mrcas[i][j][t] = mrca;

                    int firstCladeSize = mrca.getChild(0).getLeafNodeCount();
                    int secondCladeSize = mrca.getChild(1).getLeafNodeCount();

                    int possibleCombinations = firstCladeSize * secondCladeSize;
                    double treeScore = 1.0 / possibleCombinations / observedTrees.size();

                    pairwiseScoresMatrix[i][j] += treeScore;
                    pairwiseScoresPerTree[i][j][t] = treeScore;
                }
            }
        }

        // build greedy configuration

        List<Pair<Integer, Integer>> specifiedDistances = new ArrayList<>();
        Set<Integer> leavesUsed = new HashSet<>();

        while (leavesUsed.size() < n) {
            // find best next pair

            double bestScore = Double.NEGATIVE_INFINITY;
            Pair<Integer, Integer> bestPair = null;

            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    if (!leavesUsed.isEmpty() && !leavesUsed.contains(i) && !leavesUsed.contains(j)) continue;
                    if (leavesUsed.contains(i) && leavesUsed.contains(j)) continue;

                    if (bestScore < pairwiseScoresMatrix[i][j]) {
                        bestScore = pairwiseScoresMatrix[i][j];
                        bestPair = new Pair<>(i, j);
                    }
                }
            }

            specifiedDistances.add(bestPair);
            leavesUsed.add(bestPair.getFirst());
            leavesUsed.add(bestPair.getSecond());

            // adapt scores

            for (int t = 0; t < observedTrees.size(); t++) {
                Node pairMRCA = mrcas[bestPair.getFirst()][bestPair.getSecond()][t];

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (i == j) continue;

                        if (pairMRCA != mrcas[i][j][t]) continue;

                        pairwiseScoresMatrix[i][j] -= pairwiseScoresPerTree[i][j][t];
                    }
                }
            }
        }

        return CubeUtils.getConfiguration(
                specifiedDistances, tree
        );
    }
}
