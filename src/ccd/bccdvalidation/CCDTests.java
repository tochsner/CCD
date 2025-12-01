package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import beast.base.parser.NexusParser;
import ccd.model.*;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CCDTests {
    static File REFERENCE_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/subsampled");
    static File QUERY_TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/beast");
    static File PROBABILITIES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/probabilities");

    private static List<Tree> loadTrees(File treeFile) throws IOException {
        NexusParser parser = new NexusParser();
        parser.parseFile(treeFile);
        return parser.trees;
    }

    public static void calculateCCDProbabilities(
            File referenceTreeFile,
            File queryTreeFile,
            File outputProbabilitiesFile,
            AbstractCCD ccd
    ) throws IOException {
        List<Tree> referenceTrees = loadTrees(referenceTreeFile);
        List<Tree> queryTrees = loadTrees(queryTreeFile);

        for (Tree tree : referenceTrees) {
            ccd.addTree(tree);
        }
        ccd.initialize();

        Map<String, Double> cache = new HashMap<>();

        List<Double> probabilities = new ArrayList<>();
        for (Tree tree : queryTrees) {
            String prunedNewick = tree.getRoot().toNewick(true);
            double treeProbability = cache.computeIfAbsent(prunedNewick, (k) -> ccd.getProbabilityOfTree(tree));
            probabilities.add(treeProbability);
        }

        try (PrintWriter out = new PrintWriter(outputProbabilitiesFile)) {
            for (int i = 0; i < queryTrees.size(); i++) {
                out.println("STATE" + queryTrees.get(i).getID() + "," + probabilities.get(i));
            }
        }
    }

    private static File getQueryFile(File referenceFile) {
        String fileName = referenceFile.getName();

        // fileName is something like yule-10-0-3.trees,
        // where yule-10-0 is the dataset and 3 the number
        // of samples in the subsample

        // the query file is then called yule-10-0.trees
        int lastHyphenIndex = fileName.lastIndexOf("-");
        String queryFileName = fileName.substring(0, lastHyphenIndex) + ".trees";

        return new File(QUERY_TREES_DIR, queryFileName);
    }

    private static File getProbabilitiesFile(File referenceFile, String method) {
        String fileName = referenceFile.getName();

        // fileName is something like yule-10-0-3.trees,
        // where yule-10-0 is the dataset and 3 the number
        // of samples in the subsample

        // the probabilities file is then called yule-10-0-3-(method).csv
        int lastDotIndex = fileName.lastIndexOf(".");
        String refFileName = fileName.substring(0, lastDotIndex) + "-" + method + ".csv";

        return new File(PROBABILITIES_DIR, refFileName);
    }

    public static void main(String[] args) throws IOException {
        System.out.println("Load trees...");

        int numProcessedTrees = 0;
        int numLeaves = 10;

        File[] referenceFiles = REFERENCE_TREES_DIR.listFiles();
        for (File referenceFile : referenceFiles) {
            if (!referenceFile.isFile() || !referenceFile.getName().endsWith(".trees")) {
                continue;
            }

            System.out.println("Process " + referenceFile.getName() + " (" + numProcessedTrees + "/" + referenceFiles.length + ")");

            File queryFile = getQueryFile(referenceFile);

            File probabilitiesFileCCD0 = getProbabilitiesFile(referenceFile, "CCD0");
            File probabilitiesFileCCD1 = getProbabilitiesFile(referenceFile, "CCD1");
            File probabilitiesFileCCD2 = getProbabilitiesFile(referenceFile, "CCD2");

//            AbstractCCD ccd0 = new CCD0(numLeaves, false);
            AbstractCCD ccd1 = new CCD1(numLeaves, false);
//            AbstractCCD ccd2 = new CCD2(numLeaves, false);

//            calculateCCDProbabilities(referenceFile, queryFile, probabilitiesFileCCD0, ccd0);
            calculateCCDProbabilities(referenceFile, queryFile, probabilitiesFileCCD1, ccd1);
//            calculateCCDProbabilities(referenceFile, queryFile, probabilitiesFileCCD2, ccd2);

            numProcessedTrees += 1;
        }
    }
}
