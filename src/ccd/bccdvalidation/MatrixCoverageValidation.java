package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.model.matrix.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

public class MatrixCoverageValidation {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/matrix_coverage_data");

    static String DATASET_NAME = "yule-50";

    static Map<String, MatrixEstimator> ESTIMATORS = Map.of(
            "random-cube", new RandomCubeEstimator(),
            "greedy-cube", new GreedyCubeEstimator(),
//            "optimal-cube", new OptimalCubeEstimator(),
//            "mst-matrix", new MSTMatrixEstimator(),
            "greedy-matrix", new GreedyMatrixEstimator()
    );

    public static void main(String[] args) throws IOException {
        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
            if (!treeFile.getName().endsWith(".trees") || !treeFile.getName().startsWith(DATASET_NAME)) continue;

            String dataSetName = getDataSetName(treeFile);
            System.out.println("Process dataset " + dataSetName + "...");

            List<Tree> trees = Utils.loadTrees(treeFile);
            Integer numLeaves = trees.get(0).getLeafNodeCount();

            try {
                for (Map.Entry<String, MatrixEstimator> estimatorEntry : ESTIMATORS.entrySet()) {
                    System.out.println("Process estimator " + estimatorEntry.getKey() + "...");

                    MatrixCCD ccd = estimatorEntry.getValue().buildCCD(numLeaves, true);
                    for (Tree tree : trees) {
                        ccd.addTree(tree);
                    }

                    TreeMatrixConfiguration configuration = estimatorEntry.getValue().estimateConfiguration(ccd);

                    int numTreesCovered = 0;
                    for (Tree tree : trees) {
                        if (configuration.isTreeCompatible(tree)) {
                            numTreesCovered += 1;
                        }
                    }

                    try (FileWriter file = new FileWriter(getCoverageFile(dataSetName, estimatorEntry.getKey()))) {
                        file.write(String.valueOf(1.0 * numTreesCovered / trees.size()));
                    }
                }
            } catch (Exception e) {
                System.out.println(e);
            }
        }
    }

    /* -- FILE NAME HELPER FUNCTIONS -- */

    private static String getDataSetName(File treeFile) {
        return Utils.getFileNameWoExtension(treeFile);
    }

    private static File getCoverageFile(String dataSetName, String estimatorName) {
        String logFileName = dataSetName + "_" + estimatorName + "_coverage.log";
        return new File(OUT_DIR, logFileName);
    }
}
