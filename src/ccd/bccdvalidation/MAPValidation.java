package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.model.*;

import java.io.*;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class MAPValidation {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/map_data");

    static String DATASET_NAME = "yule-20";

    static int ALL_SAMPLES = -1;
    static int[] SAMPLE_SIZES = new int[]{
            10,
            100,
            1000,
            ALL_SAMPLES
    };

    static Map<String, ParameterEstimator> ESTIMATORS = Map.of(
//            "gamma-corr", new BCCDGammaMLE(new BCCDCorrelationFeatureSelector()),
//            "log-normal-corr", new BCCDLogNormalMLE(new BCCDCorrelationFeatureSelector(true)),
//            "dirichlet", new SBCCDDirichlet(),
//            "sb-beta", new SBCCDIndependentBeta(),
//            "sb-beta-per-clade", new SBCCDCladeBeta()
//            "ctmc-exp-clade", new CTMCCCDExponential(),
//            "ctmc-gamma-clade", new CTMCCCDGamma()
//            "ind-matrix-mean", new FullMatrixEstimator()
            "logit-multivariate-gaussian", new SBCCDLogisticMultivariateGaussian()
    );

    public static void main(String[] args) throws IOException {
        if (args.length == 3) {
            DATASET_NAME = args[0];
            MCMC_TREES_DIR = new File(args[1]);
            OUT_DIR = new File(args[2]);
        }

        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
            if (!treeFile.getName().startsWith(DATASET_NAME) || !treeFile.getName().endsWith(".trees")) continue;

            String dataSetName = getDataSetName(treeFile);
            System.out.println("Process dataset " + dataSetName + "...");

            List<Tree> trees = Utils.loadTrees(treeFile);
            int numLeaves = trees.get(0).getLeafNodeCount();

            // we shuffle the trees to be able to efficiently sample by taking the first n elements
            Collections.shuffle(trees);

            for (int sampleSize : SAMPLE_SIZES) {
                if (sampleSize == ALL_SAMPLES) sampleSize = trees.size();

                for (Map.Entry<String, ParameterEstimator> estimatorEntry : ESTIMATORS.entrySet()) {
                    String modelName = estimatorEntry.getKey();
                    ParameterEstimator estimator = estimatorEntry.getValue();

                    System.out.println("Run model " + modelName + " for " + sampleSize + " samples...");

                    AbstractCCD ccd = estimator.buildCCD(numLeaves, true);
                    for (int i = 0; i < sampleSize; i++) {
                        ccd.addTree(trees.get(i));
                    }
                    ccd.initialize();

                    Tree BCCDTree = ccd.getMAPTree(HeightSettingStrategy.None);

                    Utils.storeTrees(List.of(BCCDTree), getMAPTreeName(dataSetName, modelName, sampleSize));
                }

                CCD1 ccd = new CCD1(numLeaves, true);
                for (int i = 0; i < sampleSize; i++) {
                    ccd.addTree(trees.get(i));
                }
                ccd.initialize();

                Tree MRCATree = ccd.getMAPTree(HeightSettingStrategy.CommonAncestorHeights);
                Utils.storeTrees(List.of(MRCATree), getMAPTreeName(dataSetName, "mrca", sampleSize));
            }
        }
    }

    /* -- FILE NAME HELPER FUNCTIONS -- */

    private static String getDataSetName(File treeFile) {
        return Utils.getFileNameWoExtension(treeFile);
    }

    private static File getMAPTreeName(String dataSetName, String modelName, int sampleSize) {
        String mapTreeFileName = dataSetName + "_" + sampleSize + "_" + modelName + ".trees";
        return new File(OUT_DIR, mapTreeFileName);
    }
}
