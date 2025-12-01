package ccd.bccdvalidation;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class CollectMetadata {
    static File TREES_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data/thinned_mcmc_runs");

    public static void main(String[] args) throws IOException {
        if (args.length == 2) {
            TREES_DIR = new File(args[0]);
            OUT_DIR = new File(args[1]);
        }

        List<File> treeFiles = new ArrayList<>();

        for (File file : TREES_DIR.listFiles()) {
            if (file.getName().endsWith(".trees")) treeFiles.add(file);

            if (file.isDirectory()) {
                for (File subFile : file.listFiles()) {
                    if (subFile.getName().endsWith(".trees")) treeFiles.add(subFile);
                }
            }
        }

        for (File treeFile : treeFiles) {
            System.out.println("Process " + treeFile.getName());

            try {
                List<Tree> trees = Utils.loadTrees(treeFile);

                if (trees == null) {
                    System.out.println("Failed loading.");
                    continue;
                }

                int numTrees = trees.size();
                int numTaxa = trees.get(0).getLeafNodeCount();

                double leafHeight = trees.get(0).getExternalNodes().get(0).getHeight();
                boolean isContemporarySampled = true;

                for (Node leaf : trees.get(0).getExternalNodes()) {
                    if (Math.abs(leafHeight - leaf.getHeight()) > 1e-8)
                        isContemporarySampled = false;
                }

                double treeHeight = trees.get(0).getRoot().getHeight() - leafHeight;

                StringBuilder output = new StringBuilder();
                output.append(treeFile.getName()).append(",");
                output.append(numTrees).append(",");
                output.append(numTaxa).append(",");
                output.append(isContemporarySampled).append(",");
                output.append(treeHeight).append("\n");
                try (FileWriter file = new FileWriter(new File(OUT_DIR, "metadata.csv"), true)) {
                    file.write(output.toString());
                }
            } catch (Exception e) {
                System.out.println("Error: " + e);
            }
        }
    }
}
