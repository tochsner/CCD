package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.model.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class EntropyValidation {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/entropy_data");

    static String DATASET_NAME = "yule-10";

    public static void main(String[] args) throws IOException {
        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
            if (!treeFile.getName().endsWith(".trees") || !treeFile.getName().startsWith(DATASET_NAME)) continue;

            String dataSetName = getDataSetName(treeFile);
            System.out.println("Process dataset " + dataSetName + "...");

            List<Tree> trees = Utils.loadTrees(treeFile);
            Integer numLeaves = trees.get(0).getLeafNodeCount();

            CCD0 ccd = new CCD0(numLeaves, true);
            for (Tree tree : trees) {
                ccd.addTree(tree);
            }
            ccd.initialize();
            ccd.expand();

            double entropy = ccd.getEntropy();

            try (FileWriter file = new FileWriter(getEntropyFile(dataSetName))) {
                file.write(String.valueOf(entropy));
            }
        }
    }

    /* -- FILE NAME HELPER FUNCTIONS -- */

    private static String getDataSetName(File treeFile) {
        return Utils.getFileNameWoExtension(treeFile);
    }

    private static File getEntropyFile(String dataSetName) {
        String logFileName = dataSetName + "_entropy.log";
        return new File(OUT_DIR, logFileName);
    }
}
