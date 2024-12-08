package ccd.model.matrix;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeUtils;
import org.apache.commons.math3.util.Pair;

import java.util.*;

public class CubeUtils {
    static TreeMatrixConfiguration getConfiguration(int[] order, Tree tree) {
        List<Pair<Integer, Integer>> distancesSpecified = new ArrayList<>();
        for (int i = 0; i < order.length - 1; i++) {
            distancesSpecified.add(new Pair<>(order[i], order[i+1]));
        }

        List<Taxon> taxaNames = new ArrayList<>();
        for (int i = 0; i < tree.getLeafNodeCount(); i++) {
            taxaNames.add(new Taxon(String.valueOf(i)));
        }
        TaxonSet taxonSet = new TaxonSet(taxaNames);

        TreeMatrixConfiguration configuration = new TreeMatrixConfiguration(
                tree.getLeafNodeCount(),
                distancesSpecified,
                taxonSet
        );
        return configuration;
    }

    static int getPathLength(Tree tree, int node1Idx, int node2Idx) {
        if (node1Idx == node2Idx) return 0;

        Node node1 = tree.getNode(node1Idx);
        Node node2 = tree.getNode(node2Idx);
        Node mrca = TreeUtils.getCommonAncestorNode(tree, Set.of(node1.getID(), node2.getID()));

        int mrcaToNode1 = 0;
        while (node1.getParent() != mrca) {
            mrcaToNode1++;
            node1 = node1.getParent();
        }

        int mrcaToNode2 = 0;
        while (node2.getParent() != mrca) {
            mrcaToNode2++;
            node2 = node2.getParent();
        }

        return 1 + mrcaToNode1 + mrcaToNode2;
    }
}
