package ccd.model.matrix;

import beast.base.evolution.tree.Tree;
import ccd.model.HeightSettingStrategy;
import org.jgrapht.Graph;
import org.jgrapht.GraphPath;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.jgrapht.alg.tour.HeldKarpTSP;

import java.util.*;

public class OptimalCubeEstimator extends MatrixEstimator {
    @Override
    TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd) {
        Tree mapTopology = ccd.getMAPTree(HeightSettingStrategy.CommonAncestorHeights);
        TreeMatrixConfiguration configuration = this.getOptimalCube(mapTopology, ccd.getBaseTrees());
        return configuration;
    }

    TreeMatrixConfiguration getOptimalCube(Tree tree, List<Tree> observedTrees) {
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
                    int pathLength = CubeUtils.getPathLength(observedTree, i, j);
                    score += Math.log(Math.pow(2, 1 - pathLength));
                }

                graph.addEdge(String.valueOf(i), String.valueOf(j));
                graph.setEdgeWeight(graph.getEdge(String.valueOf(i), String.valueOf(j)), -score);

                pairwiseScoresMatrix[i][j] = score;
            }
        }

        GraphPath cycle = new HeldKarpTSP().getTour(graph);

        LinkedList<Integer> order = new LinkedList<>();
        for (Object vertex : cycle.getVertexList()) {
            order.add(Integer.parseInt((String) vertex));
        }
        order.remove(n);

        this.optimizeOrderToBeCompatible(tree, order, n, pairwiseScoresMatrix);

        return CubeUtils.getConfiguration(
                order.stream().mapToInt(i -> i).toArray(), tree
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
