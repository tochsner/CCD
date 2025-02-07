package ccd.model.matrix;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeUtils;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.Graph;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.interfaces.ShortestPathAlgorithm;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.SimpleWeightedGraph;
import org.jgrapht.alg.shortestpath.BFSShortestPath;

import java.util.*;

public class CubeUtils {
    static TreeMatrixConfiguration getConfiguration(int[] order, Tree tree) {
        List<Pair<Integer, Integer>> distancesSpecified = new ArrayList<>();
        for (int i = 0; i < order.length - 1; i++) {
            distancesSpecified.add(new Pair<>(order[i], order[i + 1]));
        }

        List<Taxon> taxaNames = new ArrayList<>();
        for (int i = 0; i < tree.getLeafNodeCount(); i++) {
            taxaNames.add(new Taxon(tree.getNode(i).getID()));
        }
        TaxonSet taxonSet = new TaxonSet(taxaNames);

        TreeMatrixConfiguration configuration = new TreeMatrixConfiguration(
                tree.getLeafNodeCount(),
                distancesSpecified,
                taxonSet
        );
        return configuration;
    }

    static TreeMatrixConfiguration getConfiguration(List<Pair<Integer, Integer>> distancesSpecified, Tree tree) {
        List<Taxon> taxaNames = new ArrayList<>();
        for (int i = 0; i < tree.getLeafNodeCount(); i++) {
            taxaNames.add(new Taxon(tree.getNode(i).getID()));
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

    static SimpleGraph createUnweightedGraphForTree(Tree tree) {
        SimpleGraph<Integer, DefaultEdge> graph = new SimpleGraph<>(DefaultEdge.class);
        createGraphForVertex(graph, tree.getRoot());
        return graph;
    }

    static void createGraphForVertex(SimpleGraph graph, Node vertex) {
        graph.addVertex(vertex.getNr());

        for (Node child : vertex.getChildren()) {
            createGraphForVertex(graph, child);
            graph.addEdge(vertex.getNr(), child.getNr());
        }
    }

    public static int[][] getDistanceMatrix(SimpleGraph<Integer, DefaultEdge> tree, int numLeaves) {
        int[][] distances = new int[numLeaves][numLeaves];

        BFSShortestPath bfs = new BFSShortestPath(tree);

        for (int i = 0; i < numLeaves; i++) {
            ShortestPathAlgorithm.SingleSourcePaths<Integer, DefaultEdge> vertexDistances = bfs.getPaths(i);

            for (int j = 0; j < numLeaves; j++) {
                distances[i][j] = (int) vertexDistances.getWeight(j);
            }
        }

        return distances;
    }

    public static int[][] getMrcas(SimpleGraph<Integer, DefaultEdge> tree, int numLeaves, int root) {
        BFSShortestPath bfs = new BFSShortestPath(tree);
        ShortestPathAlgorithm.SingleSourcePaths<Integer, DefaultEdge> pathsFromRoot = bfs.getPaths(root);

        int[][] mrcas = new int[numLeaves][numLeaves];

        for (int i = 0; i < numLeaves; i++) {
            List<Integer> pathFromRoot1 = pathsFromRoot.getPath(i).getVertexList();

            for (int j = i + 1; j < numLeaves; j++) {
                List<Integer> pathFromRoot2 = pathsFromRoot.getPath(j).getVertexList();

                for (int k = 0; k < pathFromRoot1.size(); k++) {
                    if (pathFromRoot1.get(k) != pathFromRoot2.get(k)) {
                        mrcas[i][j] = pathFromRoot1.get(k - 1);
                        mrcas[j][i] = pathFromRoot1.get(k - 1);
                        break;
                    }
                }
            }
        }

        return  mrcas;
    }
}
