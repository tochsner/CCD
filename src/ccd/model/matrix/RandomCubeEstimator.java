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

public class RandomCubeEstimator extends MatrixEstimator {
    static Random random = new Random();

    @Override
    TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd) {
        Tree mapTopology = ccd.getMAPTree(HeightSettingStrategy.CommonAncestorHeights);
        TreeMatrixConfiguration configuration = this.getRandomCompatibleCube(mapTopology);
        return configuration;
    }

    TreeMatrixConfiguration getRandomCompatibleCube(Tree tree) {
        int[] order = new int[tree.getLeafNodeCount()];
        int[] runningIndex = new int[]{0};

        getRandomCompatibleCube(tree.getRoot(), order, runningIndex);

        return CubeUtils.getConfiguration(order, tree);
    }

    private int[] getRandomCompatibleCube(Node vertex, int[] order, int[] runningIndex) {
        if (vertex.isLeaf()) {
            order[runningIndex[0]++] = vertex.getNr();
        } else {
            int firstIndex = random.nextInt(2);

            Node firstChild = vertex.getChild(firstIndex);
            Node secondChild = vertex.getChild(1 - firstIndex);

            getRandomCompatibleCube(firstChild, order, runningIndex);
            getRandomCompatibleCube(secondChild, order, runningIndex);
        }

        return order;
    }
}
