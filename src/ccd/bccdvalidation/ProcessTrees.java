package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.TreeDistances;
import ccd.model.*;
import ccd.tools.TraceStatistics;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.util.Map.entry;

public class ProcessTrees {
    static File MCMC_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs_replica");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/complete_data");

    static String DATASET_NAME = "yule-10";

    static Map<String, ParameterEstimator> ESTIMATORS = Map.ofEntries(
            // entry("gamma", new BCCDGammaMoM()),
            // entry("gamma-corr", new BCCDGammaMLE(new BCCDCorrelationFeatureSelector())),
            // entry("log-normal", new BCCDLogNormalMoM()),
            // entry("log-normal-corr", new BCCDLogNormalMLE(new BCCDCorrelationFeatureSelector(true))),
            // entry("sb-beta", new SBCCDIndependentBeta()),
            entry("sb-beta-per-clade", new SBCCDCladeBeta())
            // entry("sb-beta-prior", new SBCCDIndependentBetaWithPrior()),
            // entry("dirichlet", new SBCCDDirichlet()),
            // entry("cube-matrix-mean", new RandomCubeEstimator()),
            // entry("greedy-matrix-mean", new GreedyMatrixEstimator())
    );

    static int NUM_TREES_TO_SAMPLE = 10_000;

    static double TRAIN_FRACTION = 0.75;
    static double BURN_IN_FRACTION = 0.1;

    static int SEED = 0;

    public static void main(String[] args) throws IOException {
        if (args.length == 3) {
            DATASET_NAME = args[0];
            MCMC_TREES_DIR = new File(args[1]);
            OUT_DIR = new File(args[2]);
        }

//        for (File treeFile : MCMC_TREES_DIR.listFiles()) {
//            if (!treeFile.getName().endsWith(".trees")) continue;
//            if (!treeFile.getName().startsWith(DATASET_NAME)) continue;
//
//            processTree(treeFile);
//        }

        processTree(new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs_replica/yule-10-0.trees"));
    }

    public static void processTree(File treeFile) throws IOException {
        String baseFileName = treeFile.getName().split("\\.")[0];

        System.out.println("Load reference trees for " + baseFileName + "...");

        List<Tree> referenceTrees = Utils.loadTrees(treeFile);
        referenceTrees = referenceTrees.subList((int) Math.floor(referenceTrees.size() * BURN_IN_FRACTION), referenceTrees.size());

        System.out.println("Get tree ESS and subsample...");

        int treeEss = getTreeESS(referenceTrees);
        int numLeaves = referenceTrees.get(0).getLeafNodeCount();

        List<Tree> treesSubsampled = referenceTrees; // thinOutTrees(referenceTrees, treeEss);

        System.out.println("Split into train and val set for validation...");

        int numTrainTrees = (int) Math.floor(treesSubsampled.size() * TRAIN_FRACTION);

        List<Tree> treesTrain = treesSubsampled.subList(0, numTrainTrees);
        List<Tree> treesVal = treesSubsampled.subList(numTrainTrees, treesSubsampled.size());

        System.out.println("Get entropy...");

        AbstractCCD ccd1 = new CCD1(numLeaves, true);
        for (Tree tree : treesTrain) {
            ccd1.addTree(tree);
        }
        ccd1.initialize();

        double entropy = ccd1.getEntropy();

        System.out.println("Get MRCA point estimate...");

        HashMap<String, Tree> pointEstimates = new HashMap<>();
        Tree mrcaTree = ccd1.getMAPTree(HeightSettingStrategy.CommonAncestorHeights);
        pointEstimates.put("MRCA", mrcaTree);

        System.out.println("Log general stats...");

        StringBuilder statsOutput = new StringBuilder();
        statsOutput.append("distribution;metric;value").append("\n");
        statsOutput.append("-;tree_ess;").append(treeEss).append("\n");
        statsOutput.append("-;entropy;").append(entropy).append("\n");
        statsOutput.append("-;num_taxa;").append(numLeaves).append("\n");
        statsOutput.append("-;num_trees;").append(referenceTrees.size()).append("\n");

        StringBuilder credibleOutput = new StringBuilder();
        credibleOutput.append("distribution;sample_credible_set").append("\n");

        for (Map.Entry<String, ParameterEstimator> estimatorEntry : ESTIMATORS.entrySet()) {
            try {
                String modelName = estimatorEntry.getKey();
                ParameterEstimator estimator = estimatorEntry.getValue();

                System.out.println("Fit distribution " + modelName + " on train set...");

                AbstractCCD ccd = estimator.buildCCD(numLeaves, true);
                for (Tree tree : treesTrain) {
                    ccd.addTree(tree);
                }
                ccd.initialize();

                System.out.println("Get density...");

                List<Double> logDensitiesVal = new ArrayList<>();
                double finiteLogDensityVal = 0.0;
                for (Tree tree : treesVal) {
                    double logDensity = ccd.getLogProbabilityOfTree(tree);
                    logDensitiesVal.add(logDensity);
                    finiteLogDensityVal += logDensity != Double.NEGATIVE_INFINITY ? logDensity : 0.0;
                }

                System.out.println("Get point estimates...");

                Tree ccdMapTree = ccd.getMAPTree(HeightSettingStrategy.None);
                pointEstimates.put(modelName, ccdMapTree);

                System.out.println("Sample trees...");

                List<Tree> samples = new ArrayList<>();
                List<Double> logDensitiesSamples = new ArrayList<>();
                for (int i = 0; i < NUM_TREES_TO_SAMPLE; i++) {
                    Tree sample = ccd.sampleTree();
                    samples.add(sample);
                    logDensitiesSamples.add(ccd.getLogProbabilityOfTree(sample));
                }

                System.out.println("Get credible sets...");

                List<Double> credibleSets = credibleSets(logDensitiesSamples, logDensitiesVal);
                double cramerVonMises = cramerVonMises(credibleSets);

                System.out.println("Log results...");

                statsOutput.append(modelName).append(";log_data_likelihood;").append(finiteLogDensityVal).append("\n");
                statsOutput.append(modelName).append(";cramer_von_mises;").append(cramerVonMises).append("\n");

                System.out.println("Log credible sets...");

                for (Double credibleSet : credibleSets) {
                    credibleOutput.append(modelName).append(";").append(credibleSet).append("\n");
                }
            } catch (Exception e) {
                System.out.println("Error: " + e.getMessage());
            }
        }

        System.out.println("Store results...");

        try (FileWriter file = new FileWriter(new File(OUT_DIR, baseFileName + "_stats.log"))) {
            file.write(statsOutput.toString());
        }

        try (FileWriter file = new FileWriter(new File(OUT_DIR, baseFileName + "_credible_sets.log"))) {
            file.write(credibleOutput.toString());
        }

        System.out.println("Store point estimates trees...");

        Utils.storeTrees(pointEstimates, new File(OUT_DIR, baseFileName + "_point_estimate.trees"));
    }

    static int getTreeESS(List<Tree> trees) throws IOException {
        CCD0 ccd = new CCD0(trees.get(0).getLeafNodeCount(), false);

        for (Tree tree : trees) {
            ccd.addTree(tree);
        }

        ccd.initialize();
        WrappedBeastTree ccdMAP = new WrappedBeastTree(ccd.getMAPTree());

        // calculate expected RF distance to CCD1 MAP

        double[] expectedRfTrace = new double[trees.size()];

        for (int i = 0; i < trees.size(); i++) {
            Tree tree = trees.get(i);
            WrappedBeastTree wrappedTree = new WrappedBeastTree(tree);
            expectedRfTrace[i] = TreeDistances.robinsonsFouldDistance(
                    wrappedTree,
                    ccdMAP
            );
        }

        // estimate ESS

        TraceStatistics traceStatistics = new TraceStatistics(expectedRfTrace, 1);
        int ess = (int) Math.floor(traceStatistics.getESS());

        return ess;
    }

    static List<Tree> thinOutTrees(List<Tree> trees, int ess) throws IOException {
        // choose which indices to sample

        List<Integer> chosenIndices = IntStream.range(0, trees.size()).boxed().collect(Collectors.toList());
        Collections.shuffle(chosenIndices, new Random(SEED));
        chosenIndices = chosenIndices.subList(0, ess);
        Collections.sort(chosenIndices);

        // sub-sample trees

        List<Tree> chosenTrees = new ArrayList<>();
        for (int i : chosenIndices) {
            chosenTrees.add(trees.get(i));
        }

        return chosenTrees;
    }

    public static List<Double> credibleSets(List<Double> references, List<Double> queries) {
        Collections.sort(references);

        List credibleSets = new ArrayList();

        for (Double query: queries) {
            int ecdf = references.size();

            for (int i = 0; i < references.size(); i++) {
                if (references.get(i) > query) {
                    ecdf = i;
                    break;
                }
            }

            credibleSets.add(1.0 - 1.0 * ecdf / references.size());
        }

        return credibleSets;
    }

    public static double cramerVonMises(List<Double> x) {
        int n = x.size();
        Collections.sort(x);
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            double term = x.get(i) - ((2.0 * (i + 1) - 1.0) / (2.0 * n));
            sum += term * term;
        }
        return (1.0 / (12.0 * n) + sum) / n;
    }

}
