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

    public SBCCD(int numLeaves, boolean storeBaseTrees, ParameterEstimator<SBCCD> estimator) {
        super(numLeaves, storeBaseTrees);
        this.estimator = estimator;
    }

    /* -- INITIALIZATION -- */

    ParameterEstimator<SBCCD> estimator;

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

    private Clade getClade(Node vertex) {
        BitSet cladeInBits = BitSet.newBitSet(leafArraySize);
        Clade firstChildClade;
        Clade secondChildClade;

        if (vertex.isLeaf()) {
            int index = vertex.getNr();
            cladeInBits.set(index);
        } else {
            firstChildClade = this.getClade(vertex.getChildren().get(0));
            secondChildClade = this.getClade(vertex.getChildren().get(1));

            cladeInBits.or(firstChildClade.getCladeInBits());
            cladeInBits.or(secondChildClade.getCladeInBits());
        }

        return this.cladeMapping.get(cladeInBits);
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
        if (heightStrategy != HeightSettingStrategy.None) {
            return super.getTreeBasedOnStrategy(samplingStrategy, heightStrategy);
        }

        double treeHeight = this.sampleTreeHeight();

        Node root = this.getVertexBasedOnStrategy(
                samplingStrategy,
                this.rootClade,
                new int[]{this.getSizeOfLeavesArray()},
                treeHeight
        );
        return new Tree(root);
    }

    protected Node getVertexBasedOnStrategy(
            SamplingStrategy samplingStrategy,
            Clade clade,
            int[] runningInnerIndex,
            double vertexHeight
    ) {
        if (clade.isLeaf()) {
            int leafNr = clade.getCladeInBits().nextSetBit(0);
            String taxonName = this.getSomeBaseTree().getTaxaNames()[leafNr];

            assert vertexHeight == 0.0;

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
            firstVertexHeight = 0.0;
        } else if (samplingStrategy == SamplingStrategy.Sampling) {
            firstVertexHeight = vertexHeight - partition.sampleFirstBranchLength(vertexHeight);
        } else if (samplingStrategy == SamplingStrategy.MAP) {
            firstVertexHeight = vertexHeight - partition.getMAPFirstBranchLength(vertexHeight);
        } else {
            throw new UnsupportedOperationException("This sampling strategy is not yet implemented.");
        }

        if (secondClade.isLeaf()) {
            secondVertexHeight = 0.0;
        } else if (samplingStrategy == SamplingStrategy.Sampling) {
            secondVertexHeight = vertexHeight - partition.sampleSecondBranchLength(vertexHeight);
        } else if (samplingStrategy == SamplingStrategy.MAP) {
            secondVertexHeight = vertexHeight - partition.getMAPSecondBranchLength(vertexHeight);
        } else {
            throw new UnsupportedOperationException("This sampling strategy is not yet implemented.");
        }

        Node firstChild = this.getVertexBasedOnStrategy(
                samplingStrategy,
                firstClade,
                runningInnerIndex,
                firstVertexHeight
        );
        Node secondChild = this.getVertexBasedOnStrategy(
                samplingStrategy,
                secondClade,
                runningInnerIndex,
                secondVertexHeight
        );

        Node vertex = new Node();
        vertex.setNr(runningInnerIndex[0]++);

        vertex.setHeight(vertexHeight);
        vertex.addChild(firstChild);
        vertex.addChild(secondChild);

        return vertex;
    }

    @Override
    public void sampleBranchLengths(Tree tree) {
        double treeHeight = this.sampleTreeHeight();

        Node root = tree.getRoot();
        root.setHeight(treeHeight);

        this.sampleBranchLengths(root);
    }

    protected void sampleBranchLengths(Node vertex) {
        double vertexHeight = vertex.getHeight();

        if (vertex.isLeaf()) {
            assert vertexHeight == 0.0;

            int index = vertex.getNr();

            BitSet cladeInBits = BitSet.newBitSet(leafArraySize);
            cladeInBits.set(index);

            return;
        }

        Node firstChild = vertex.getChild(0);
        Node secondChild = vertex.getChild(1);

        Clade clade = this.getClade(vertex);
        Clade firstClade = this.getClade(firstChild);
        Clade secondClade = this.getClade(secondChild);

        SBCCDCladePartition partition = (SBCCDCladePartition) clade.getCladePartition(firstClade);

        double firstVertexHeight;
        double secondVertexHeight;

        if (firstClade.isLeaf()) {
            firstVertexHeight = 0.0;
        } else {
            firstVertexHeight = vertexHeight - partition.getMAPFirstBranchLength(vertexHeight);
        }

        if (secondClade.isLeaf()) {
            secondVertexHeight = 0.0;
        } else {
            secondVertexHeight = vertexHeight - partition.getMAPSecondBranchLength(vertexHeight);
        }

        if (firstVertexHeight < 0) {
            throw new RuntimeException();
        }

        firstChild.setHeight(firstVertexHeight);
        secondChild.setHeight(secondVertexHeight);

        sampleBranchLengths(firstChild);

        sampleBranchLengths(secondChild);
    }

    /* -- TREE LIKELIHOOD - TREE LIKELIHOOD -- */

    @Override
    public double getProbabilityOfTree(Tree tree) {
        double height = SBCCDCladePartition.getMaxDistanceToLeaf(tree.getRoot());
        double heightDensity = this.getHeightDistribution().density(height);

        double treeDensity = super.getProbabilityOfTree(tree);

        return heightDensity * treeDensity;
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
