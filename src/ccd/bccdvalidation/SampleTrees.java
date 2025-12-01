package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.model.AbstractCCD;
import ccd.model.CCD1;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class SampleTrees {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    static File TRUE_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_config");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/ccd1_sample_data");

    static String DATASET_NAME = "yule-50_94";

    static int NUM_TREES_TO_PROCESS = 100;
    static int NUM_TREES_TO_SAMPLE = 10_000;

    public static void main(String[] args) throws IOException {
        int numTreesProcessed = 0;

        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
            if (!treeFile.getName().startsWith(DATASET_NAME) || !treeFile.getName().endsWith(".trees")) continue;

            String dataSetName = Utils.getFileNameWoExtension(treeFile);
            System.out.println("Process dataset " + dataSetName + "...");

            List<Tree> trees = Utils.loadTrees(treeFile);
            int numLeaves = trees.get(0).getLeafNodeCount();

            // build CCD based on the MCMC trees

            AbstractCCD ccd = new CCD1(numLeaves, true);
            for (Tree tree : trees) {
                ccd.addTree(tree);
            }
            ccd.initialize();

            // calculate the posterior of the true tree

            StringBuilder output = new StringBuilder();
            output.append("tree,log_posterior\n");

            Tree trueTree = Utils.loadTrees(getTrueTreeFile(treeFile)).get(0);
            double trueTreePosterior = ccd.getLogProbabilityOfTree(trueTree);
            output.append("true,");
            output.append(trueTreePosterior);
            output.append("\n");

            // calculate the posterior of trees sampled with the BCCD model

            List<Tree> sampledTrees = new ArrayList<>();
            for (int i = 0; i < NUM_TREES_TO_SAMPLE; i++) {
                Tree sampledTree = ccd.sampleTree();
                double sampleTreePosterior = ccd.getLogProbabilityOfTree(sampledTree);

                output.append(i);
                output.append(",");
                output.append(sampleTreePosterior);
                output.append("\n");

                sampledTrees.add(sampledTree);
            }

            try (FileWriter file = new FileWriter(getPosteriorFileName(dataSetName))) {
                file.write(output.toString());
            }

            // store trees

            sampledTrees.add(0, trueTree);
            Utils.storeTrees(sampledTrees, getTreesFileName(dataSetName));

            if (++numTreesProcessed == NUM_TREES_TO_PROCESS) break;
        }
    }

    /* -- FILE NAME HELPER FUNCTIONS -- */

    private static File getTrueTreeFile(File treeFile) {
        return new File(TRUE_TREES_DIR, treeFile.getName());
    }

    private static File getPosteriorFileName(String dataSetName) {
        String mapTreeFileName = dataSetName + ".log";
        return new File(OUT_DIR, mapTreeFileName);
    }
    private static File getTreesFileName(String dataSetName) {
        String mapTreeFileName = dataSetName + ".trees";
        return new File(OUT_DIR, mapTreeFileName);
    }
}
