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
import org.jgrapht.graph.SimpleWeightedGraph;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

public class MSTMatrixEstimator extends MatrixEstimator {
    @Override
    public TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd) {
        Tree mapTopology = ccd.getMAPTree(HeightSettingStrategy.CommonAncestorHeights);
        TreeMatrixConfiguration configuration = this.getOptimalConfiguration(mapTopology, ccd.getBaseTrees());
        return configuration;
    }

    TreeMatrixConfiguration getOptimalConfiguration(Tree tree, List<Tree> observedTrees) {
        int n = tree.getLeafNodeCount();

        Graph<String, DefaultEdge> graph = new SimpleWeightedGraph<>(DefaultEdge.class);
        for (int i = 0; i < n; i++) {
            graph.addVertex(String.valueOf(i));
        }

        double[][] pairwiseScoresMatrix = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) continue;

                double score = 0;

                for (Tree observedTree : observedTrees) {
                    Node node1 = observedTree.getNode(i);
                    Node node2 = observedTree.getNode(j);
                    Node mrca = TreeUtils.getCommonAncestorNode(
                            observedTree,
                            Set.of(node1.getID(), node2.getID())
                    );

                    int firstCladeSize = mrca.getChild(0).getLeafNodeCount();
                    int secondCladeSize = mrca.getChild(1).getLeafNodeCount();

                    int possibleCombinations = firstCladeSize * secondCladeSize;
                    score += Math.log(1.0 / possibleCombinations);
                }

                graph.addEdge(String.valueOf(i), String.valueOf(j));
                graph.setEdgeWeight(graph.getEdge(String.valueOf(i), String.valueOf(j)), -score);

                pairwiseScoresMatrix[i][j] = score;
            }
        }

        SpanningTreeAlgorithm.SpanningTree cycle = new PrimMinimumSpanningTree(graph).getSpanningTree();

        List<Pair<Integer, Integer>> specifiedDistances = new ArrayList<>();
        for (Object edge : cycle.getEdges()) {
            specifiedDistances.add(new Pair<>(
                    Integer.parseInt(graph.getEdgeSource((DefaultEdge) edge)),
                    Integer.parseInt(graph.getEdgeTarget((DefaultEdge) edge))
            ));
        }

        return CubeUtils.getConfiguration(
                specifiedDistances, tree
        );
    }

    private static void optimizeOrderToBeCompatible(Tree tree, LinkedList<Integer> greedyOrder, int n, double[][] pairwiseScoresMatrix) {
        int bestFirstLeaf = greedyOrder.getFirst();
        double bestRotationScore = Double.NEGATIVE_INFINITY;

        for (int offset = 0; offset < n; offset++) {
            TreeMatrixConfiguration configuration = CubeUtils.getConfiguration(
                    greedyOrder.stream().mapToInt(i -> i).toArray(), tree
            );
            if (!configuration.isTreeCompatible(tree)) {
                // rotate by one
                greedyOrder.addLast(greedyOrder.pollFirst());
                continue;
            }

            double offsetScore = 0;
            for (int i = 0; i < n - 1; i++) {
                int firstLeaf = greedyOrder.get(i);
                int secondLeaf = greedyOrder.get(i + 1);
                offsetScore += pairwiseScoresMatrix[firstLeaf][secondLeaf];
            }

            if (bestRotationScore < offsetScore) {
                bestRotationScore = offsetScore;
                bestFirstLeaf = greedyOrder.getFirst();
            }

            // rotate by one
            greedyOrder.addLast(greedyOrder.pollFirst());
        }

        while (greedyOrder.getFirst() != bestFirstLeaf) {
            greedyOrder.addLast(greedyOrder.pollFirst());
        }
    }
}
