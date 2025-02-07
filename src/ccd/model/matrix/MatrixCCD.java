package ccd.model.matrix;

import beast.base.evolution.tree.Tree;
import ccd.model.*;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

public class MatrixCCD extends AbstractCCD {
    TreeMatrixConfiguration configuration;
    MultivariateNormalDistribution cubeDistribution;
    ParameterEstimator estimator;

    public MatrixCCD(int numLeaves, boolean storeTrees, ParameterEstimator estimator) {
        super(numLeaves, storeTrees);
        this.estimator = estimator;
    }

    @Override
    public void initialize() {
        this.estimator.estimateParameters(this);
    }

    public TreeMatrixConfiguration getConfiguration() {
        return configuration;
    }

    /** -- PARAMETERS - PARAMETERS -- */

    public void setConfiguration(TreeMatrixConfiguration configuration) {
        this.configuration = configuration;
    }

    public MultivariateNormalDistribution getCubeDistribution() {
        return cubeDistribution;
    }

    public void setCubeDistribution(MultivariateNormalDistribution cubeDistribution) {
        this.cubeDistribution = cubeDistribution;
    }

    /* -- PROBABILITY - PROBABILITY -- */

    @Override
    public double getProbabilityOfTree(Tree tree) {
        throw new UnsupportedOperationException();
    }

    public double getLogProbabilityOfTree(Tree tree) {
        if (!this.getConfiguration().isTreeCompatible(tree)) return Double.NEGATIVE_INFINITY;

        double[] cube = this.getConfiguration().getDistancesForCompatibleTree(tree);

        double logCubeProbability = Math.log(this.getCubeDistribution().density(cube));
        return logCubeProbability;
    }

    /* -- TREE SAMPLING - TREE SAMPLING -- */

    @Override
    public Tree sampleTree() {
        double[] cube = this.getCubeDistribution().sample();
        return this.getConfiguration().getTree(cube);
    }

    /* -- MAP TREE -- */

    @Override
    public Tree getMAPTree(HeightSettingStrategy settingStrategy) {
        if (settingStrategy != HeightSettingStrategy.None) return super.getMAPTree(settingStrategy);

        double[] cube = this.getCubeDistribution().getMeans();
        return this.getConfiguration().getTree(cube);
    }

    /** -- UNSUPPORTED METHODS **/

    @Override
    protected boolean removeCladePartitionIfNecessary(Clade clade, CladePartition partition) {
        throw new UnsupportedOperationException();
    }

    @Override
    protected void tidyUpCacheIfDirty() {
        // nothing to do
    }

    @Override
    public double getNumberOfParameters() {
        throw new UnsupportedOperationException();
    }

    @Override
    public AbstractCCD copy() {
        throw new UnsupportedOperationException();
    }
}
