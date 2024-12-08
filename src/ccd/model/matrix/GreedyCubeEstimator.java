package ccd.model.matrix;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.HeightSettingStrategy;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class GreedyCubeEstimator extends MatrixEstimator {
    @Override
    TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd) {
        Tree mapTopology = ccd.getMAPTree(HeightSettingStrategy.CommonAncestorHeights);
        TreeMatrixConfiguration configuration = this.getGreedyCube(mapTopology, ccd.getBaseTrees());
        return configuration;
    }

    TreeMatrixConfiguration getGreedyCube(Tree tree, List<Tree> observedTrees) {
        int n = tree.getLeafNodeCount();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double score = 0;

                for (Tree observedTree : observedTrees) {
                    int pathLength = CubeUtils.getPathLength(observedTree, i, j);
                    score +=
                }
            }
        }

        return null;
    }
}
