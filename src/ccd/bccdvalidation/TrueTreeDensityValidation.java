package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import ccd.model.*;
//import org.apache.commons.math3.stat.descriptive.rank.Median;
//import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.io.*;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class TrueTreeDensityValidation {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs");
    static File TRUE_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_config");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/true_tree_density_data");

    static String DATASET_NAME = "yule-50";

    static Map<String, ParameterEstimator> ESTIMATORS = Map.of(
//            "gamma-corr", new BCCDGammaMLE(new BCCDCorrelationFeatureSelector()),
//            "log-normal", new BCCDLogNormalMoM()
//            "log-normal-corr", new BCCDLogNormalMLE(new BCCDCorrelationFeatureSelector(true)),
            "dirichlet-sample", new SBCCDDirichlet(18.2405410265742)
//            "sb-beta-per-clade-sample", new SBCCDCladeBeta()
//            "ccd0", new CCD1Estimator()
    );

    static int NUM_TREES_TO_PROCESS = 10000;
    static int NUM_TREES_TO_SAMPLE = 10_000;

    public static void main(String[] args) throws IOException {
        System.out.println("Run true tree density validation with " + args.length + " args.");

        if (args.length == 4) {
            DATASET_NAME = args[0];
            MCMC_TREES_DIR = new File(args[1]);
            TRUE_TREES_DIR = new File(args[2]);
            OUT_DIR = new File(args[3]);
        }

        System.out.println("Start");

        long s = System.nanoTime();
        for (int i = 0; i < 20; i++) {
                TreeAnnotator.FastTreeSet treeSet = new TreeAnnotator().new FastTreeSet("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs/yule-50_98.trees", 0);
        }
        long e = System.nanoTime();
        System.out.println((e - s) / 1_000_000);

        int numTreesProcessed = 0;

        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
            if (!treeFile.getName().startsWith(DATASET_NAME) || !treeFile.getName().endsWith(".trees")) continue;

            String dataSetName = Utils.getFileNameWoExtension(treeFile);
            System.out.println("Process dataset " + dataSetName + "...");

            List<Tree> trees = Utils.loadTrees(new File("/Users/tobiaochsner/Documents/Thesis/TractableTreeDistributions/test/ref_trees.trees"));
//            List<Tree> trees = ccd.bccdvalidation.Utils.loadTrees(treeFile);
            List<Tree> query_trees = Utils.loadTrees(new File("/Users/tobiaochsner/Documents/Thesis/TractableTreeDistributions/test/query_trees.trees"));
            int numLeaves = trees.get(0).getLeafNodeCount();

            Random random = new Random();

            Tree trueTree = Utils.loadTrees(getTrueTreeFile(treeFile)).get(0);

            for (Map.Entry<String, ParameterEstimator> estimatorEntry : ESTIMATORS.entrySet()) {
                String modelName = estimatorEntry.getKey();
                ParameterEstimator estimator = estimatorEntry.getValue();

                File outputFile = getPosteriorFileName(dataSetName, modelName);

                System.out.println("Run model " + modelName + "...");

//                try {
                // build BCCD based on the MCMC trees

                long start = System.nanoTime();

                AbstractCCD ccd = estimator.buildCCD(numLeaves, true);
                for (Tree tree : trees) {
                    ccd.addTree(tree);
                }
                ccd.initialize();

                long end = System.nanoTime();
                System.out.println((end - start) / 1_000_000);
                ccd.getLogProbabilityOfTree(query_trees.get(0));
                // calculate the posterior of the true tree

                StringBuilder output = new StringBuilder();
                output.append("tree,log_posterior\n");

                double trueTreePosterior = ccd.getLogProbabilityOfTree(trees.get(random.nextInt(trees.size())));
                output.append("true,");
                output.append(trueTreePosterior);
                output.append("\n");

                // calculate the posterior of trees sampled with the BCCD model

                for (int i = 0; i < NUM_TREES_TO_SAMPLE; i++) {
                    Tree sampledTree = ccd.sampleTree();
                    double sampleTreePosterior = ccd.getLogProbabilityOfTree(sampledTree);

                    output.append(i);
                    output.append(",");
                    output.append(sampleTreePosterior);
                    output.append("\n");
                }

                    try (FileWriter file = new FileWriter(outputFile)) {
                        file.write(output.toString());
                    }

//                } catch (Exception exp) {
//                    System.out.println("Error: " + exp);
//                }
            }

            if (++numTreesProcessed == NUM_TREES_TO_PROCESS) break;
        }
    }

    /* -- FILE NAME HELPER FUNCTIONS -- */

    private static File getTrueTreeFile(File treeFile) {
        return new File(TRUE_TREES_DIR, treeFile.getName());
    }

    private static File getPosteriorFileName(String dataSetName, String modelName) {
        // dataSetName is something like yule-10_0,
        // the output file is then called yule-10_0_<modelName>.log
        String mapTreeFileName = dataSetName + "_" + modelName + ".log";
        return new File(OUT_DIR, mapTreeFileName);
    }

}
