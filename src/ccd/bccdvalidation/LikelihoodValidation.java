package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.model.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class LikelihoodValidation {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/likelihood_data");

    static String DATASET_NAME = "yule-20";
    static int ALL_SAMPLES = -1;
    static int[] SAMPLE_SIZES = new int[] {
//            10,
//            100,
//            1000,
            ALL_SAMPLES
    };
    static int SEED = 100;

    static Map<String, ParameterEstimator> ESTIMATORS = Map.of(
            "gamma-beta-corr", new BCCDGammaMLE(new BCCDCovarianceFeatureSelector(false)),
            "log-normal-beta-corr", new BCCDLogNormalMLE(new BCCDCorrelationFeatureSelector(true)),
            "sb-dirichlet", new SBCCDDirichlet(),
            "sb-beta", new SBCCDCladeBeta()
    );

    public static void main(String[] args) throws IOException {
        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
            if (!treeFile.getName().endsWith(".trees") || !treeFile.getName().startsWith(DATASET_NAME)) continue;

            String dataSetName = getDataSetName(treeFile);
            System.out.println("Process dataset " + dataSetName + "...");

            List<Tree> trees = Utils.loadTrees(treeFile);
            Integer numLeaves = trees.get(0).getLeafNodeCount();

            // we shuffle the trees such that we can easily obtain random samples by
            // taking the first n elements
            Collections.shuffle(trees, new Random(SEED));

            for (int sampleSize : SAMPLE_SIZES) {
                if (sampleSize == ALL_SAMPLES) {
                    sampleSize = trees.size();
                }

                for (Map.Entry<String, ParameterEstimator> estimatorEntry : ESTIMATORS.entrySet()) {
                    String modelName = estimatorEntry.getKey();
                    System.out.println("Run model " + modelName + "...");

                    ParameterEstimator estimator = estimatorEntry.getValue();
                    AbstractCCD ccd = estimator.buildCCD(numLeaves, false);

                    for (int i = 0; i < sampleSize; i++) {
                        ccd.addTree(trees.get(i));
                    }
                    ccd.initialize();

                    calculateMCMCTreePosteriors(ccd, trees, dataSetName, modelName, sampleSize);
                }
            }
        }
    }

    /* -- TESTS -- */

    private static void calculateMCMCTreePosteriors(AbstractCCD ccd, List<Tree> referenceTrees, String dataSetName, String modelName, int sampleSize) throws IOException {
        StringBuilder output = new StringBuilder();

        output.append("state,log_posterior,num_parameters\n");

        for (Tree tree : referenceTrees) {
            double posterior = ccd.getLogProbabilityOfTree(tree);
            output.append("STATE" + tree.getID() + ",");
            output.append(posterior);
            output.append(",");
            output.append(ccd.getNumberOfParameters());
            output.append("\n");
        }

        try (FileWriter file = new FileWriter(getPosteriorsFile(dataSetName, modelName, sampleSize))) {
            file.write(output.toString());
        }
    }

    /* -- FILE NAME HELPER FUNCTIONS -- */

    private static String getDataSetName(File treeFile) {
        return Utils.getFileNameWoExtension(treeFile);
    }

    private static File getPosteriorsFile(String dataSetName, String modelName, int sampleSize) {
        String logFileName = dataSetName + "_" + sampleSize + "_logs_" + modelName + ".log";
        return new File(OUT_DIR, logFileName);
    }
}
