package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import beast.base.parser.NexusParser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;

public class Utils {
    public static String getFileNameWoExtension(File file) {
        String fileName = file.getName();
        int lastPeriodIndex = fileName.lastIndexOf(".");
        return fileName.substring(0, lastPeriodIndex);
    }

    public static List<Tree> loadTrees(File treeFile) throws IOException {
        NexusParser parser = new NexusParser();
        parser.parseFile(treeFile);
        return parser.trees;
    }

    public static void storeTrees(List<Tree> trees, File treeFile) throws FileNotFoundException {
        try (PrintStream handle = new PrintStream(treeFile)) {
            trees.get(0).init(handle);

            for (int i = 0; i < trees.size(); i++) {
                Tree tree = trees.get(i);

                int id = (tree.getID() == null) ? i : Integer.parseInt(tree.getID().replace("_", ""));
                tree.log(id, handle);
                handle.println();
            }

            trees.get(trees.size() - 1).close(handle);
        }
    }

    public static void storeTrees(Map<String, Tree> trees, File treeFile) throws FileNotFoundException {
        try (PrintStream handle = new PrintStream(treeFile)) {
            trees.values().stream().findFirst().orElseThrow().init(handle);
            handle.println();

            for (Map.Entry<String, Tree> entry : trees.entrySet()) {
                String label = entry.getKey();
                Tree tree = entry.getValue();

                handle.print("tree '" + label + "' = ");

                final String newick = tree.getRoot().toSortedNewick(new int[1], false);
                handle.print(newick);
                handle.print(";");

                handle.println();
            }

            trees.values().stream().findFirst().orElseThrow().close(handle);
        }
    }
}
