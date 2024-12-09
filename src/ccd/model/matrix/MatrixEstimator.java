package ccd.model.matrix;

import beast.base.evolution.tree.Tree;
import ccd.model.HeightSettingStrategy;
import ccd.model.ParameterEstimator;
import ccd.tools.MatrixUtils;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;

import java.util.ArrayList;
import java.util.List;

public abstract class MatrixEstimator extends ParameterEstimator<MatrixCCD> {
    @Override
    public MatrixCCD buildCCD(int numLeaves, boolean storeBaseTrees) {
        if (!storeBaseTrees) throw new IllegalArgumentException("Base trees must be stored to use the CubeEstimator");
        return new MatrixCCD(numLeaves, storeBaseTrees, this);
    }

    @Override
    public void estimateParameters(MatrixCCD ccd) {
        TreeMatrixConfiguration configuration = this.estimateConfiguration(ccd);
        MultivariateNormalDistribution cubeDistribution = this.estimateDistribution(ccd, configuration);

        ccd.setConfiguration(configuration);
        ccd.setCubeDistribution(cubeDistribution);
    }

    public abstract TreeMatrixConfiguration estimateConfiguration(MatrixCCD ccd);

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

        double[] means = new double[ccd.getNumberOfLeaves() - 1];
        for (double[] cube : cubes) {
            for (int i = 0; i < cube.length; i++) {
                means[i] += cube[i] / cubes.length;
            }
        }

        RealMatrix covarianceMatrix = new Covariance(cubes).getCovarianceMatrix();
        double[][] covariances = new double[covarianceMatrix.getRowDimension()][covarianceMatrix.getColumnDimension()];
        MatrixUtils.fillArray(covarianceMatrix, covariances);

        MultivariateNormalDistribution distribution = new MultivariateNormalDistribution(means, covariances);
        return distribution;
    }

    @Override
    public int getNumberOfParameters(MatrixCCD ccd) {
        throw new UnsupportedOperationException();
    }
}
