package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.model.*;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class DistributionValidation {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/distribution_data");

    static String DATASET_NAME = "yule-50_94";

    static Map<String, ParameterEstimator> ESTIMATORS = Map.of(
//            "log-normal-mu-sigma", new BCCDLogNormalMoM(),
//            "log-normal-beta-old-old", new BCCDLogNormalMLE(CladePartitionObservation::branchLengthOldOld),
//            "log-normal-beta-cov", new BCCDLogNormalMLE(new BCCDCovarianceFeatureSelector(true)),
//            "log-normal-beta-corr", new BCCDLogNormalMLE(new BCCDCorrelationFeatureSelector(true)),
//            "gamma-mu-sigma", new BCCDGammaMoM(),
//            "gamma-beta-old-old", new BCCDGammaMLE(CladePartitionObservation::branchLengthOldOld),
//            "gamma-beta-cov", new BCCDGammaMLE(new BCCDCovarianceFeatureSelector()),
            "gamma-beta-corr", new BCCDGammaMLE(new BCCDCorrelationFeatureSelector()),
//            "dirichlet", new SBCCDDirichlet()
            "sb-beta", new SBCCDIndependentBeta()
    );

    static int NUM_TREES_TO_SAMPLE = 5_000;

    public static void main(String[] args) throws IOException {
        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
            if (!treeFile.getName().endsWith(".trees") || !treeFile.getName().startsWith(DATASET_NAME)) continue;

            String dataSetName = getDataSetName(treeFile);
            System.out.println("Process dataset " + dataSetName + "...");

            List<Tree> trees = Utils.loadTrees(treeFile);
            Integer numLeaves = trees.get(0).getLeafNodeCount();

            for (Map.Entry<String, ParameterEstimator> estimatorEntry : ESTIMATORS.entrySet()) {
                String modelName = estimatorEntry.getKey();
                ParameterEstimator estimator = estimatorEntry.getValue();

                System.out.println("Run model " + modelName + "...");

                AbstractCCD ccd = estimator.buildCCD(numLeaves, false);
                for (Tree tree : trees) {
                    ccd.addTree(tree);
                }
                ccd.initialize();

                sampleTrees(ccd, dataSetName, modelName);
                calculateMCMCTreePosteriors(ccd, trees, dataSetName, modelName);
            }
        }
    }

    /* -- TESTS -- */

    private static void sampleTrees(
            AbstractCCD ccd, String dataSetName, String modelName
    ) throws FileNotFoundException {
        List<Tree> sampledTrees = new ArrayList<>();

        for (int i = 0; i < NUM_TREES_TO_SAMPLE; i++) {
            Tree sampledTree = ccd.sampleTree(HeightSettingStrategy.None);
            sampledTrees.add(sampledTree);
        }

        Utils.storeTrees(sampledTrees, getSampledTreesFile(dataSetName, modelName));
    }

    private static void calculateMCMCTreePosteriors(AbstractCCD ccd, List<Tree> referenceTrees, String dataSetName, String modelName) throws IOException {
        StringBuilder output = new StringBuilder();

        output.append("state,posterior\n");

        for (Tree tree : referenceTrees) {
            double posterior = ccd.getProbabilityOfTree(tree);
            output.append("STATE" + tree.getID() + ",");
            output.append(posterior);
            output.append("\n");
        }

        try (FileWriter file = new FileWriter(getPosteriorsFile(dataSetName, modelName))) {
            file.write(output.toString());
        }
    }

    /* -- FILE NAME HELPER FUNCTIONS -- */

    private static String getDataSetName(File treeFile) {
        return Utils.getFileNameWoExtension(treeFile);
    }

    private static File getSampledTreesFile(String dataSetName, String modelName) {
        // dataSetName is something like yule-10-0,
        // the sampled trees file is then called yule-10-0_sampled-trees_<modelName>.trees
        String sampledTreesFileName = dataSetName + "_sampled-trees_" + modelName + ".trees";
        return new File(OUT_DIR, sampledTreesFileName);
    }

    private static File getSampledBranchLengthsFile(String dataSetName, String modelName) {
        // fileName is something like yule-10-0.trees,
        // the sampled trees file is then called yule-10-0_sampled-branches_<modelName>.trees
        String sampledTreesFileName = dataSetName + "_sampled-branches_" + modelName + ".trees";
        return new File(OUT_DIR, sampledTreesFileName);
    }

    private static File getPosteriorsFile(String dataSetName, String modelName) {
        // fileName is something like yule-10-0.trees,
        // the log file is then called yule-10-0-<modelName>.log
        String logFileName = dataSetName + "_logs_" + modelName + ".log";
        return new File(OUT_DIR, logFileName);
    }
}
