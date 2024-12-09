package ccd.model.matrix;

import beast.base.evolution.tree.Tree;
import ccd.model.HeightSettingStrategy;
import org.apache.commons.math3.util.Pair;

import java.util.*;

public class GreedyCubeEstimator extends MatrixEstimator {
    @Override
    TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd) {
        Tree mapTopology = ccd.getMAPTree(HeightSettingStrategy.CommonAncestorHeights);
        TreeMatrixConfiguration configuration = this.getGreedyCube(mapTopology, ccd.getBaseTrees());
        return configuration;
    }

    TreeMatrixConfiguration getGreedyCube(Tree tree, List<Tree> observedTrees) {
        int n = tree.getLeafNodeCount();

        double[][] pairwiseScoresMatrix = new double[n][n];
        SortedMap<Double, Integer>[] sortedPairwiseScores = new SortedMap[n];

        double bestScore = Double.NEGATIVE_INFINITY;
        Pair<Integer, Integer> bestPair = new Pair<>(0, 1);

        for (int i = 0; i < n; i++) {
            sortedPairwiseScores[i] = new TreeMap<>(Collections.reverseOrder());

            for (int j = 0; j < n; j++) {
                if (i == j) continue;

                double score = 0;

                for (Tree observedTree : observedTrees) {
                    int pathLength = CubeUtils.getPathLength(observedTree, i, j);
                    score += Math.log(Math.pow(2, 1 - pathLength));
                }

                if (bestScore < score) {
                    bestScore = score;
                    bestPair = new Pair<>(i, j);
                }

                pairwiseScoresMatrix[i][j] = score;

                // we add a small jitter to make the keys in the sorted map unique
                score *= 1 + 1e-6 * j;
                sortedPairwiseScores[i].put(score, j);
            }
        }

        LinkedList<Integer> greedyOrder = getGreedyOrder(bestPair, n, sortedPairwiseScores);
        optimizeGreedyOrderToBeCompatible(tree, greedyOrder, n, pairwiseScoresMatrix);

        return CubeUtils.getConfiguration(
                greedyOrder.stream().mapToInt(i -> i).toArray(), tree
        );
    }

    private LinkedList<Integer> getGreedyOrder(Pair<Integer, Integer> bestPair, int n, SortedMap<Double, Integer>[] sortedPairwiseScores) {
        LinkedList<Integer> greedyOrder = new LinkedList<>();

        greedyOrder.add(bestPair.getFirst());
        greedyOrder.add(bestPair.getSecond());

        while (greedyOrder.size() < n) {
            int currentHead = greedyOrder.getFirst();
            int currentTail = greedyOrder.getLast();

            removeExistingScores(
                    sortedPairwiseScores[currentHead], greedyOrder
            );
            removeExistingScores(
                    sortedPairwiseScores[currentTail], greedyOrder
            );

            double bestNewHeadScore = sortedPairwiseScores[currentHead].firstKey();
            double bestNewTailScore = sortedPairwiseScores[currentTail].firstKey();

            if (bestNewTailScore < bestNewHeadScore) {
                greedyOrder.addFirst(sortedPairwiseScores[currentHead].get(bestNewHeadScore));
            } else {
                greedyOrder.addLast(sortedPairwiseScores[currentTail].get(bestNewTailScore));
            }
        }
        return greedyOrder;
    }

    private void removeExistingScores(SortedMap<Double, Integer> sortedScores, List<Integer> existingOrder) {
        Set<Integer> leavesInExistingOrder = new HashSet<>();
        for (int leaf : existingOrder) {
            leavesInExistingOrder.add(leaf);
        }

        while (!sortedScores.isEmpty()) {
            double bestScore = sortedScores.firstKey();
            int bestLeaf = sortedScores.get(bestScore);

            if (!leavesInExistingOrder.contains(bestLeaf)) {
                return;
            }

            sortedScores.remove(bestScore);
        }
    }

    private static void optimizeGreedyOrderToBeCompatible(Tree tree, LinkedList<Integer> greedyOrder, int n, double[][] pairwiseScoresMatrix) {
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
