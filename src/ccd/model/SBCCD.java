package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

public class SBCCD extends AbstractCCD {

    /* -- CONSTRUCTORS & CONSTRUCTION METHODS -- */

    public SBCCD(int numLeaves, boolean storeBaseTrees, SBCCDParameterEstimator estimator) {
        super(numLeaves, storeBaseTrees);
        this.estimator = estimator;
    }

    /* -- INITIALIZATION -- */

    SBCCDParameterEstimator estimator;

    @Override
    public void initialize() {
        this.estimator.estimateParameters(this);
    }

    @Override
    protected void initializeRootClade(int numLeaves) {
        this.leafArraySize = numLeaves;

        BitSet rootBitSet = BitSet.newBitSet(leafArraySize);
        rootBitSet.set(0, numLeaves);

        this.rootClade = new SBCCDClade(rootBitSet, this);
        cladeMapping.put(rootClade.getCladeInBits(), rootClade);
    }

    /* -- CLADE CREATION - CLADE CREATION -- */

    @Override
    protected Clade addNewClade(BitSet cladeInBits) {
        Clade clade = new SBCCDClade(cladeInBits, this);
        cladeMapping.put(cladeInBits, clade);
        return clade;
    }

    /* -- HEIGHT SAMPLING - HEIGHT SAMPLING -- */

    protected BranchLengthDistribution heightDistribution;

    public BranchLengthDistribution getHeightDistribution() {
        return heightDistribution;
    }

    public void setHeightDistribution(BranchLengthDistribution heightDistribution) {
        this.heightDistribution = heightDistribution;
    }

    protected double sampleTreeHeight() {
        BranchLengthDistribution dist = this.getHeightDistribution();
        return dist.sample();
    }

    /* -- TREE SAMPLING - TREE SAMPLING -- */

    @Override
    protected Tree getTreeBasedOnStrategy(SamplingStrategy samplingStrategy, HeightSettingStrategy heightStrategy) {
        if (samplingStrategy != SamplingStrategy.Sampling || heightStrategy != HeightSettingStrategy.None) {
            return super.getTreeBasedOnStrategy(samplingStrategy, heightStrategy);
        }

        Node root = this.getVertexBasedOnStrategy(
                this.rootClade,
                new int[]{this.getSizeOfLeavesArray()},
                0,
                this.sampleTreeHeight()
        );
        return new Tree(root);
    }

    protected Node getVertexBasedOnStrategy(
            Clade clade,
            int[] runningInnerIndex,
            double vertexHeight,
            double subTreeHeight
    ) {
        if (clade.isLeaf()) {
            int leafNr = clade.getCladeInBits().nextSetBit(0);
            String taxonName = this.getSomeBaseTree().getTaxaNames()[leafNr];

            assert subTreeHeight == 0.0;

            Node vertex = new Node(taxonName);
            vertex.setNr(leafNr);
            vertex.setHeight(vertexHeight);

            return vertex;
        }

        SBCCDCladePartition partition = (SBCCDCladePartition) this.getPartitionBasedOnStrategy(clade, SamplingStrategy.Sampling);
        if (partition == null) {
            throw new AssertionError("Unsuccessful to find clade partition of clade: " + clade.getCladeInBits());
        }

        Clade firstClade = partition.getChildClades()[0];
        Clade secondClade = partition.getChildClades()[1];

        double firstVertexHeight;
        double secondVertexHeight;

        if (firstClade.isLeaf()) {
            firstVertexHeight = vertexHeight + subTreeHeight;
        } else {
            firstVertexHeight = partition.sampleFirstBranchLength(subTreeHeight);
        }

        if (secondClade.isLeaf()) {
            secondVertexHeight = vertexHeight + subTreeHeight;
        } else {
            secondVertexHeight = partition.sampleSecondBranchLength(subTreeHeight);
        }

        Node firstChild = this.getVertexBasedOnStrategy(
                firstClade,
                runningInnerIndex,
                firstVertexHeight,
                subTreeHeight - firstVertexHeight
        );
        Node secondChild = this.getVertexBasedOnStrategy(
                secondClade,
                runningInnerIndex,
                secondVertexHeight,
                subTreeHeight - secondVertexHeight
        );

        Node vertex = new Node();
        vertex.setNr(runningInnerIndex[0]++);

        vertex.setHeight(vertexHeight);
        vertex.addChild(firstChild);
        vertex.addChild(secondChild);

        return vertex;
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

    public List<SBCCDCladePartition> getAllPartitions() {
        Set<SBCCDCladePartition> partitions = new HashSet<>();

        for (Clade clade : this.getClades()) {
            for (CladePartition partition : clade.partitions) {
                partitions.add((SBCCDCladePartition) partition);
            }
        }

        return new ArrayList<>(partitions);
    }

    @Override
    public String toString() {
        return "SBCCD " + super.toString();
    }

    @Override
    public AbstractCCD copy() {
       throw new UnsupportedOperationException();
    }

    @Override
    protected double getNumberOfParameters() {
        return 3 * this.getNumberOfCladePartitions();
    }

}
