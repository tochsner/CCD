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

        Node root = this.getVertexBasedOnStrategy(
                samplingStrategy,
                this.rootClade,
                new int[]{this.getSizeOfLeavesArray()},
                0.0
        );

        double timeSinceLastSplit = this.sampleTimeSinceLastSplit(samplingStrategy);
        double mostRecentSampleTime = this.getMostRecentSampleTime(root);
        this.incrementHeight(root, timeSinceLastSplit - mostRecentSampleTime);
        this.setLeavesToZero(root);

        return new Tree(root);
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

    protected Node getVertexBasedOnStrategy(
            SamplingStrategy samplingStrategy,
            Clade clade,
            int[] runningInnerIndex,
            double parentHeight
    ) {
        if (clade.isLeaf()) {
            int leafNr = clade.getCladeInBits().nextSetBit(0);
            String taxonName = this.getSomeBaseTree().getTaxaNames()[leafNr];

            Node vertex = new Node(taxonName);
            vertex.setNr(leafNr);
            vertex.setHeight(parentHeight);

            return vertex;
        }

        CTMCCCDCladePartition chosenPartition;
        double chosenRate = 0.0;

        List<CTMCCCDCladePartition> cladePartitions = clade.getPartitions().stream().map(x -> (CTMCCCDCladePartition) x).toList();
        switch (samplingStrategy) {
            case Sampling -> {
                chosenPartition = cladePartitions.get(0);

                double totalRate = cladePartitions.stream().mapToDouble(x -> x.getRate()).sum();
                double sampleAtAccumulativeRate = totalRate * random.nextDouble();

                double accumulativeRate = 0.0;
                for (CTMCCCDCladePartition partition : cladePartitions) {
                    accumulativeRate += partition.getRate();

                    if (sampleAtAccumulativeRate < accumulativeRate) {
                        chosenPartition = partition;
                        break;
                    }
                }

                chosenRate = new ExponentialDistribution(totalRate).sample();
            }
            case MAP -> {
                chosenPartition = cladePartitions.get(0);

                for (CTMCCCDCladePartition partition : cladePartitions) {
                    if (chosenRate < partition.getRate()) {
                        chosenPartition = partition;
                    }
                    chosenRate += partition.getRate();
                }
            }
            default -> {
                throw new UnsupportedOperationException("Sampling type not implemented.");
            }
        }

        double branchLength = 1 / chosenRate;

        Node firstChild = this.getVertexBasedOnStrategy(
                samplingStrategy,
                chosenPartition.getChildClades()[0],
                runningInnerIndex,
                parentHeight - branchLength
        );
        Node secondChild = this.getVertexBasedOnStrategy(
                samplingStrategy,
                chosenPartition.getChildClades()[1],
                runningInnerIndex,
                parentHeight - branchLength
        );

        Node vertex = new Node();
        vertex.setNr(runningInnerIndex[0]++);

        vertex.setHeight(parentHeight - branchLength);
        vertex.addChild(firstChild);
        vertex.addChild(secondChild);

        return vertex;
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
