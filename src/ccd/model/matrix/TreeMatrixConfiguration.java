package ccd.model.matrix;

import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.tree.ClusterTree;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeUtils;
import org.apache.commons.math3.util.Pair;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

public class TreeMatrixConfiguration {
    int n;
    TaxonSet taxonSet;
    List<Pair<Integer, Integer>> distancesSpecified;

    public TreeMatrixConfiguration(int n, List<Pair<Integer, Integer>> distancesSpecified, TaxonSet taxonSet) {
        this.n = n;
        this.distancesSpecified = distancesSpecified;
        this.taxonSet = taxonSet;
    }

    double[] getDistancesForCompatibleTree(Tree tree) {
        double[] distances = new double[n - 1];

        for (int i = 0; i < n - 1; i++) {
            Pair<Integer, Integer> distancePair = this.distancesSpecified.get(i);
            Node node1 = tree.getNode(distancePair.getFirst());
            Node node2 = tree.getNode(distancePair.getSecond());

            Node mrca = TreeUtils.getCommonAncestorNode(tree, Set.of(node1.getID(), node2.getID()));
            double logDistance = Math.log(2 * mrca.getHeight());
            distances[i] = logDistance;
        }

        return distances;
    }

    Tree getTree(double[] distances) {
        double [] distanceMatrix = new double[n * n];
        Arrays.fill(distanceMatrix, Double.MAX_VALUE);

        for (int i = 0; i < n - 1; i++) {
            Pair<Integer, Integer> distancePair = this.distancesSpecified.get(i);
            int a = distancePair.getFirst();
            int b = distancePair.getSecond();

            distanceMatrix[a * n + b] = Math.exp(distances[i]);
            distanceMatrix[a + b * n] = Math.exp(distances[i]);
        }

        Distance distance = (taxon1, taxon2) -> distanceMatrix[taxon1 * n + taxon2];

        ClusterTree clusterTree = new ClusterTree();
        clusterTree.initByName(
                "clusterType", "single",
                "distance", distance,
                "taxonset", this.taxonSet
        );

        return clusterTree;
    }

    /**
     * Return true if the tree is compatible with the tree matrix configuration.
     * <br>
     * A tree is compatible if all pairs of leaves with a specified distance
     * have a unique MRCA.
     */
    boolean isTreeCompatible(Tree tree) {
        boolean[] usedAsMRCA = new boolean[tree.getNodeCount()];

        for (Pair<Integer, Integer> distancePair : this.distancesSpecified) {
            Node node1 = tree.getNode(distancePair.getFirst());
            Node node2 = tree.getNode(distancePair.getSecond());

            Node mrca = TreeUtils.getCommonAncestorNode(tree, Set.of(node1.getID(), node2.getID()));
            int mrcaNr = mrca.getNr();

            if (usedAsMRCA[mrcaNr]) {
                return false;
            } else {
                usedAsMRCA[mrcaNr] = true;
            }
        }

        return true;
    }
}
