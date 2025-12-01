package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.algorithms.TreeDistances;
import ccd.model.CCD0;
import ccd.model.WrappedBeastTree;
import ccd.tools.TraceStatistics;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class ThinOutMCMCRuns {
    public static File COMPLETE_RUNS_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/mcmc_runs");
    public static File THINNED_RUNS_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    public static int SEED = 200;

    public static void main(String[] args) throws IOException {
        if (args.length == 2) {
            COMPLETE_RUNS_DIR = new File(args[0]);
            THINNED_RUNS_DIR = new File(args[1]);
        }

        int progress = 0;
        for (File file : COMPLETE_RUNS_DIR.listFiles()) {
            if (file.getName().endsWith(".trees")) {
                System.out.println(++progress + ": thin out " + file.getName());
                thinOutMCMCRun(file);
            }
        }
    }

    static void thinOutMCMCRun(File treesFile) throws IOException {
        if ((new File(THINNED_RUNS_DIR, treesFile.getName())).exists()) return;

        // calculate CCD0 based on trees

        List<Tree> trees = Utils.loadTrees(treesFile);
        if (trees == null) return;

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

        // choose which indices to sample

        List<Integer> chosenIndices = IntStream.range(0, trees.size()).boxed().collect(Collectors.toList());
        Collections.shuffle(chosenIndices, new Random(SEED));
        chosenIndices = chosenIndices.subList(0, ess);
        Collections.sort(chosenIndices);

        // sub-sample .trees files

        List<Tree> chosenTrees = new ArrayList<>();
        for (int i : chosenIndices) {
            chosenTrees.add(trees.get(i));
        }
        Utils.storeTrees(chosenTrees, new File(THINNED_RUNS_DIR, treesFile.getName()));

        // sub-sample .log files

        File logFile = getLogsFile(treesFile);

        List<String> relevantLines;
        try (Stream<String> lines = Files.lines(logFile.toPath())) {
            relevantLines = lines.filter(x -> !x.startsWith("#")).toList();
        }
        assert relevantLines.size() == trees.size() + 1; // first line is the header

        File thinnedLogFile = new File(THINNED_RUNS_DIR, logFile.getName());
        try (PrintStream handle = new PrintStream(thinnedLogFile)) {
            handle.println(relevantLines.get(0));   // first line is the header
            for (int i : chosenIndices) {
                handle.println(relevantLines.get(i + 1));
            }
        }
    }

    /* -- FILE NAME HELPER FUNCTIONS -- */

    private static File getLogsFile(File treeFile) {
        String treeFileName = treeFile.getName();
        int lastPeriodIndex = treeFileName.lastIndexOf(".");
        String logFileName = treeFileName.substring(0, lastPeriodIndex) + ".log";
        return new File(
                treeFile.getParent(),
                logFileName
        );
    }
}
