package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.model.*;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class DistancesExperiment {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    static File TRUE_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_config");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/distances_data");

    static String DATASET_NAME = "yule-";

    public static void main(String[] args) throws IOException {
        if (args.length == 4) {
            DATASET_NAME = args[0];
            MCMC_TREES_DIR = new File(args[1]);
            TRUE_TREES_DIR = new File(args[2]);
            OUT_DIR = new File(args[3]);
        }

        StringBuilder output = new StringBuilder();

        output.append("tree").append(",");
        output.append("avg_pairwise_distance").append(",");
        output.append("avg_min_pairwise_distance").append(",");
        output.append("avg_max_pairwise_distance").append(",");
        output.append("avg_var_pairwise_distance").append("\n");

        try (FileWriter file = new FileWriter(new File(OUT_DIR, "rf_distances.csv"))) {
            file.write(output.toString());
        }

        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
            if (!treeFile.getName().endsWith(".trees") || !treeFile.getName().startsWith(DATASET_NAME)) continue;

            System.out.println("Process " + treeFile.getName());

            System.out.println("Loading trees...");

            List<Tree> trees = Utils.loadTrees(treeFile);
            Collections.shuffle(trees);
            trees = trees.subList(0, 500);

            System.out.println("Get pair-wise distances...");

            int n = trees.size();

            // get pair-wise distances

            List<ArrayList<BitSet>> nonTrivialClades = new ArrayList<>();

            for (int i = 0; i < n; i++) {
                nonTrivialClades.add(
                    new WrappedBeastTree(trees.get(i)).getNontrivialClades()
                );
            }

            double[][] distances = new double[n][n];
            int numPairs = n * (n-1) / 2;

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    ArrayList<BitSet> firstClades = (ArrayList<BitSet>) nonTrivialClades.get(i).clone();
                    ArrayList<BitSet> secondClades = nonTrivialClades.get(j);
                    firstClades.removeAll(secondClades);
                    distances[i][j] = firstClades.size();
                }
            }

            System.out.println("Get aggregates...");

            // get average pair-wise distance

            double totalPairwiseDistance = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    totalPairwiseDistance += distances[i][j];
                }
            }
            double avgPairwiseDistance = totalPairwiseDistance / numPairs;

            // get average variance in pair-wise distances

            double totalVarPairwiseDistance = 0.0;
            for (int i = 0; i < n; i++) {
                totalVarPairwiseDistance += (new Variance()).evaluate(distances[i]);
            }
            double avgVarPairwiseDistance = totalVarPairwiseDistance / numPairs;

            // get average minimum pair-wise distance

            double totalMinPairwiseDistance = 0.0;
            for (int i = 0; i < n; i++) {
                double minDistance = Double.MAX_VALUE;
                for (int j = 0; j < n; j++) {
                    if (i != j)
                        minDistance = Math.min(minDistance, distances[i][j]);
                }
                totalMinPairwiseDistance += minDistance;
            }
            double avgMinPairwiseDistance = totalMinPairwiseDistance / n;

            // get average max pair-wise distance

            double totalMaxPairwiseDistance = 0.0;
            for (int i = 0; i < n; i++) {
                double maxDistance = Double.MIN_VALUE;
                for (int j = 0; j < n; j++) {
                    if (i != j)
                        maxDistance = Math.max(maxDistance, distances[i][j]);
                }
                totalMaxPairwiseDistance += maxDistance;
            }
            double avgMaxPairwiseDistance = totalMaxPairwiseDistance / n;

            // write output

            output = new StringBuilder();

            output.append(treeFile.getName()).append(",");
            output.append(avgPairwiseDistance).append(",");
            output.append(avgMinPairwiseDistance).append(",");
            output.append(avgMaxPairwiseDistance).append(",");
            output.append(avgVarPairwiseDistance).append("\n");

            try (FileWriter file = new FileWriter(new File(OUT_DIR, "rf_distances.csv"), true)) {
                file.write(output.toString());
            }
        }
    }

    private static File getTrueTreeFile(File treeFile) {
        return new File(TRUE_TREES_DIR, treeFile.getName());
    }
}
