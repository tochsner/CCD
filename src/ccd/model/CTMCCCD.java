package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import org.apache.commons.math3.distribution.ExponentialDistribution;

import java.util.*;

public class CTMCCCD extends AbstractCCD {

    /* -- CONSTRUCTORS & CONSTRUCTION METHODS -- */

    public CTMCCCD(int numLeaves, boolean storeBaseTrees, ParameterEstimator<CTMCCCD> estimator) {
        super(numLeaves, storeBaseTrees);
        this.estimator = estimator;
    }

    /* -- INITIALIZATION -- */

    ParameterEstimator<CTMCCCD> estimator;

    @Override
    public void initialize() {
        this.estimator.estimateParameters(this);
    }

    @Override
    protected void initializeRootClade(int numLeaves) {
        this.leafArraySize = numLeaves;

        BitSet rootBitSet = BitSet.newBitSet(leafArraySize);
        rootBitSet.set(0, numLeaves);

        this.rootClade = new CTMCCCDClade(rootBitSet, this);
        cladeMapping.put(rootClade.getCladeInBits(), rootClade);
    }

    /* -- CLADE CREATION - CLADE CREATION -- */

    @Override
    protected Clade addNewClade(BitSet cladeInBits) {
        Clade clade = new CTMCCCDClade(cladeInBits, this);
        cladeMapping.put(cladeInBits, clade);
        return clade;
    }

    /* -- HEIGHT SAMPLING - HEIGHT SAMPLING -- */

    protected BranchLengthDistribution timeSinceLastSplitDistribution;

    public BranchLengthDistribution getTimeSinceLastSplitDistribution() {
        return timeSinceLastSplitDistribution;
    }

    public void setTimeSinceLastSplitDistribution(BranchLengthDistribution timeSinceLastSplitDistribution) {
        this.timeSinceLastSplitDistribution = timeSinceLastSplitDistribution;
    }

    protected double sampleTimeSinceLastSplit(SamplingStrategy samplingStrategy) {
        return switch (samplingStrategy) {
            case Sampling -> this.getTimeSinceLastSplitDistribution().sample();
            case MAP -> {
                try {
                    yield this.getTimeSinceLastSplitDistribution().mode();
                } catch (NoModeException e) {
                    yield this.getTimeSinceLastSplitDistribution().mean();
                }
            }
            default -> throw new UnsupportedOperationException();
        };
    }

    /* -- TREE SAMPLING - TREE SAMPLING -- */

    @Override
    protected Tree getTreeBasedOnStrategy(SamplingStrategy samplingStrategy, HeightSettingStrategy heightStrategy) {
        if (heightStrategy != HeightSettingStrategy.None) {
            throw new UnsupportedOperationException();
        }

        Tree topology = super.getTreeBasedOnStrategy(samplingStrategy, heightStrategy);
        this.assignHeights(topology, samplingStrategy);

        return topology;
    }

    private void assignHeights(Tree tree, SamplingStrategy samplingStrategy) {
        Node root = tree.getRoot();

        Map<Node, CTMCCCDCladePartition> partitionsInTree = new HashMap<>();
        this.collectCladePartitions(root, partitionsInTree);

        Map<Node, Double> topBranchLengths = switch (samplingStrategy) {
            case MAP -> this.estimator.getAllMAPLengths(partitionsInTree);
            case Sampling -> this.estimator.sampleAllLengths(partitionsInTree);
            default -> throw new UnsupportedOperationException("This sampling type is not implemented.");
        };

        this.assignHeights(root, topBranchLengths, 0.0);

        double timeSinceLastSplit = this.sampleTimeSinceLastSplit(samplingStrategy);
        double mostRecentSampleTime = this.getMostRecentSampleTime(root);
        this.incrementHeight(root, timeSinceLastSplit - mostRecentSampleTime);
        this.setLeavesToZero(root);
    }

    public Clade collectCladePartitions(Node vertex, Map<Node, CTMCCCDCladePartition> collectedPartitions) {
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize);

        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            cladeInBits.set(index);
            return cladeMapping.get(cladeInBits);
        }

        Clade firstChildClade = collectCladePartitions(vertex.getChildren().get(0), collectedPartitions);
        Clade secondChildClade = collectCladePartitions(vertex.getChildren().get(1), collectedPartitions);

        cladeInBits.or(firstChildClade.getCladeInBits());
        cladeInBits.or(secondChildClade.getCladeInBits());

        Clade currentClade = cladeMapping.get(cladeInBits);

        CTMCCCDCladePartition partition = (CTMCCCDCladePartition) currentClade.getCladePartition(firstChildClade, secondChildClade);
        collectedPartitions.put(vertex, partition);

        return currentClade;
    }

    public Clade collectClades(Node vertex, Map<Node, CTMCCCDClade> collectedClades) {
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize);

        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            cladeInBits.set(index);
            return cladeMapping.get(cladeInBits);
        }

        Clade firstChildClade = collectClades(vertex.getChildren().get(0), collectedClades);
        Clade secondChildClade = collectClades(vertex.getChildren().get(1), collectedClades);

        cladeInBits.or(firstChildClade.getCladeInBits());
        cladeInBits.or(secondChildClade.getCladeInBits());

        CTMCCCDClade currentClade = (CTMCCCDClade) cladeMapping.get(cladeInBits);
        collectedClades.put(vertex, currentClade);

        return currentClade;
    }

    private void setLeavesToZero(Node vertex) {
        if (vertex.isLeaf()) {
            vertex.setHeight(0.0);
        } else {
            for (Node child : vertex.getChildren()) {
                this.setLeavesToZero(child);
            }
        }
    }

    private double getMostRecentSampleTime(Node vertex) {
        return vertex.getChildren().stream().mapToDouble(x -> getMostRecentSampleTime(x)).min().orElse(vertex.getHeight());
    }

    private void incrementHeight(Node vertex, double time) {
        vertex.setHeight(vertex.getHeight() + time);

        for (Node child : vertex.getChildren()) {
            this.incrementHeight(child, time);
        }
    }

    protected void assignHeights(
            Node vertex,
            Map<Node, Double> topBranchLengths,
            double parentHeight
    ) {
        if (vertex.isLeaf()) {
            vertex.setHeight(parentHeight);
            return;
        }

        double topBranchLength = topBranchLengths.get(vertex);
        vertex.setHeight(parentHeight - topBranchLength);

        for (Node child : vertex.getChildren()) {
            assignHeights(child, topBranchLengths, parentHeight - topBranchLength);
        }
    }

    @Override
    public void sampleBranchLengths(Tree tree) {
        throw new UnsupportedOperationException();
    }

    /* -- TREE LIKELIHOOD - TREE LIKELIHOOD -- */

    @Override
    public double getLogProbabilityOfTree(Tree tree) {
        double timeSinceLastSplit = Utils.getTimeSinceLastSplit(tree.getRoot());
        double logProbabilityTimeSinceLastSplit = this.getTimeSinceLastSplitDistribution().logDensity(
                timeSinceLastSplit
        );

        double logTreeDensity = super.getLogProbabilityOfTree(tree);

        return logProbabilityTimeSinceLastSplit + logTreeDensity;
    }

    /* -- STATE MANAGEMENT - STATE MANAGEMENT -- */

    @Override
    public void setCacheAsDirty() {
        super.setCacheAsDirty();
    }

    @Override
    protected void tidyUpCacheIfDirty() {
        resetCacheIfProbabilitiesDirty();
    }

    @Override
    protected boolean removeCladePartitionIfNecessary(Clade clade, CladePartition partition) {
        throw new UnsupportedOperationException();
    }

    /* -- OTHER METHODS -- */

    public List<CTMCCCDCladePartition> getAllPartitions() {
        Set<CTMCCCDCladePartition> partitions = new HashSet<>();

        for (Clade clade : this.getClades()) {
            for (CladePartition partition : clade.partitions) {
                partitions.add((CTMCCCDCladePartition) partition);
            }
        }

        return new ArrayList<>(partitions);
    }

    @Override
    public AbstractCCD copy() {
        throw new UnsupportedOperationException();
    }

    @Override
    public double getNumberOfParameters() {
        return this.estimator.getNumberOfParameters(this);
    }

}
