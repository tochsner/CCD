package ccd.bccdvalidation;

import beast.base.evolution.tree.Tree;
import ccd.model.AbstractCCD;
import ccd.model.CCD0;
import ccd.model.CCD1;
import ccd.model.CCD2;
import com.phylodata.loader.ExperimentLoader;
import com.phylodata.types.PaperWithExperiment;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class RankUniformityTest {
    static File MCMC_TREE_FILE = new File("/Users/tobiaochsner/Downloads/stree.trees");
    static File OUT_DIR = new File("/Users/tobiaochsner/Documents/Thesis/Validation/data");

    static int NUM_TREES_TO_SAMPLE = 10_000;

    public static void main(String[] args) throws IOException {
        PaperWithExperiment experiment = new ExperimentLoader(
                "munro-2019-climate-6tvf", 1
        ).load();

        System.out.println("Title: " + experiment.getPaper().getTitle());

//
//
//        File outputFile = new File(OUT_DIR, MCMC_TREE_FILE.getName() + ".log");
//
//        List<Tree> mcmcTrees = ccd.bccdvalidation.Utils.loadTrees(MCMC_TREE_FILE);
//
//        long start = System.nanoTime();
//        AbstractCCD ccd = new CCD1(mcmcTrees, 0.0);
//        long end = System.nanoTime();
//        System.out.println((end - start) / 1_000_000_000.0);
//
//        StringBuilder output = new StringBuilder();
//        output.append("type,probability\n");
//
//        // calculate the probability of the MCMC trees
//
//        for (Tree tree : mcmcTrees) {
//            double probability = ccd.getProbabilityOfTree(tree);
//
//            output.append("mcmc").append(",");
//            output.append(probability).append("\n");
//        }
//
//        // calculate the probability of trees sampled under CCD
//
//        for (int i = 0; i < NUM_TREES_TO_SAMPLE; i++) {
//            Tree sampledTree = ccd.sampleTree();
//            double probability = ccd.getProbabilityOfTree(sampledTree);
//
//            output.append("sampled").append(",");
//            output.append(probability).append("\n");
//        }
//
//        try (FileWriter file = new FileWriter(outputFile)) {
//            file.write(output.toString());
//        }
    }

}
