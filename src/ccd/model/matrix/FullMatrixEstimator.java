package ccd.model.matrix;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import ccd.model.HeightSettingStrategy;
import ccd.tools.MatrixUtils;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class FullMatrixEstimator extends MatrixEstimator {
    @Override
    public TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd) {
        List<Pair<Integer, Integer>> allDistances = new ArrayList<>();

        for (int i = 0; i < ccd.getNumberOfLeaves(); i++) {
            for (int j = i + 1; j < ccd.getNumberOfLeaves(); j++) {
                allDistances.add(new Pair<>(i, j));
            }
        }

        TreeMatrixConfiguration configuration = CubeUtils.getConfiguration(
                allDistances,
                ccd.getMAPTree()
        );
        configuration.allTreesAreCompatible = true;
        return configuration;
    }

    MultivariateNormalDistribution estimateDistribution(MatrixCCD ccd, TreeMatrixConfiguration configuration) {
        List<Tree> compatibleTrees = new ArrayList<>();
        for (Tree tree : ccd.getBaseTrees()) {
            if (configuration.isTreeCompatible(tree)) {
                compatibleTrees.add(tree);
            }
        }

        double[][] cubes = new double[compatibleTrees.size()][];
        for (int i = 0; i < compatibleTrees.size(); i++) {
            Tree tree = compatibleTrees.get(i);
            cubes[i] = configuration.getDistancesForCompatibleTree(tree);
        }

        double[] means = new double[configuration.distancesSpecified.size()];
        for (double[] cube : cubes) {
            for (int i = 0; i < cube.length; i++) {
                means[i] += cube[i] / cubes.length;
            }
        }

        // we assume no covariances
        double[][] covariances = new double[configuration.distancesSpecified.size()][configuration.distancesSpecified.size()];

        Variance varianceEstimator = new Variance();
        for (int i = 0; i < configuration.distancesSpecified.size(); i++) {
            covariances[i][i] = varianceEstimator.evaluate(cubes[i]);
        }

        MultivariateNormalDistribution distribution = new MultivariateNormalDistribution(means, covariances);
        return distribution;
    }
}
