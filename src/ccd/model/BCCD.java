package ccd.model;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class BCCD extends AbstractCCD {

    /* -- CONSTRUCTORS & CONSTRUCTION METHODS -- */

    public BCCD(int numLeaves, boolean storeBaseTrees, BCCDParameterEstimator estimator) {
        super(numLeaves, storeBaseTrees);
        this.estimator = estimator;
    }

    /* -- INITIALIZATION -- */

    BCCDParameterEstimator estimator;

    @Override
    public void initialize() {
        this.estimator.estimateParameters(this.getAllPartitions());
    }

    @Override
    protected void initializeRootClade(int numLeaves) {
        this.leafArraySize = numLeaves;

        BitSet rootBitSet = BitSet.newBitSet(leafArraySize);
        rootBitSet.set(0, numLeaves);

        this.rootClade = new BCCDClade(rootBitSet, this);
        cladeMapping.put(rootClade.getCladeInBits(), rootClade);
    }

    /* -- CLADE CREATION - CLADE CREATION -- */

    @Override
    protected Clade addNewClade(BitSet cladeInBits) {
        Clade clade = new BCCDClade(cladeInBits, this);
        cladeMapping.put(cladeInBits, clade);
        return clade;
    }

    /* -- TREE SAMPLING - TREE SAMPLING -- */

    @Override
    public Tree sampleTree() {
        Tree topology = sampleTree(HeightSettingStrategy.None);
        return topology;
    }

    protected Node getVertexBasedOnStrategy(
            Clade clade,
            SamplingStrategy samplingStrategy,
            HeightSettingStrategy heightStrategy
    ) {
        return this.getVertexBasedOnStrategy(
                clade,
                samplingStrategy,
                heightStrategy,
                new int[]{this.getSizeOfLeavesArray()}
        );
    }

    protected Node getVertexBasedOnStrategy(
            Clade clade,
            SamplingStrategy samplingStrategy,
            HeightSettingStrategy heightStrategy,
            int[] runningInnerIndex
    ) {
        if (heightStrategy != HeightSettingStrategy.None) {
            return super.getVertexBasedOnStrategy(clade, samplingStrategy, heightStrategy);
        }

        if (clade.isLeaf()) {
            int leafNr = clade.getCladeInBits().nextSetBit(0);
            String taxonName = this.getSomeBaseTree().getTaxaNames()[leafNr];

            Node vertex = new Node(taxonName);
            vertex.setNr(leafNr);
            vertex.setHeight(clade.getMeanOccurredHeight());

            return vertex;
        }

        BCCDCladePartition partition = (BCCDCladePartition) this.getPartitionBasedOnStrategy(clade, samplingStrategy);
        if (partition == null) {
            throw new AssertionError("Unsuccessful to find clade partition of clade: " + clade.getCladeInBits());
        }

        Node firstChild = getVertexBasedOnStrategy(partition.getChildClades()[0],
                samplingStrategy, heightStrategy);
        Node secondChild = getVertexBasedOnStrategy(partition.getChildClades()[1],
                samplingStrategy, heightStrategy);

        Node vertex = new Node();
        vertex.setNr(runningInnerIndex[0]++);

        vertex.addChild(firstChild);
        vertex.addChild(secondChild);

        double firstChildHeight = firstChild.getHeight();
        double secondChildHeight = secondChild.getHeight();

        double minBranchLength;
        if (samplingStrategy == SamplingStrategy.Sampling) {
            minBranchLength = partition.sampleMinBranchLength(vertex);
        } else if (samplingStrategy == SamplingStrategy.MAP) {
            minBranchLength = partition.getMAPBranchLength(vertex);
        } else {
            throw new AssertionError("This sampling strategy is not supported by BCCD.");
        }

        double vertexHeight = Math.max(
            firstChildHeight + minBranchLength,
            secondChildHeight + minBranchLength
        );
        vertex.setHeight(vertexHeight);

        return vertex;
    }

    public void sampleBranchLengths(Tree tree) {
        this.sampleBranchLengths(tree.getRoot());
    }

    protected Clade sampleBranchLengths(Node vertex) {
        if (vertex.isLeaf()) {
            int index = vertex.getNr();

            BitSet cladeInBits = BitSet.newBitSet(leafArraySize);
            cladeInBits.set(index);

            Clade clade = this.cladeMapping.get(cladeInBits);
            return clade;
        }

        Node firstChild = vertex.getChild(0);
        Node secondChild = vertex.getChild(1);

        Clade firstChildClade = this.sampleBranchLengths(firstChild);
        Clade secondChildClade = this.sampleBranchLengths(secondChild);

        BitSet currentCladeInBits = BitSet.newBitSet(leafArraySize);
        currentCladeInBits.or(firstChildClade.getCladeInBits());
        currentCladeInBits.or(secondChildClade.getCladeInBits());

        Clade currentClade = cladeMapping.get(currentCladeInBits);
        BCCDCladePartition currentPartition = (BCCDCladePartition) currentClade.getCladePartition(firstChildClade);

        double firstChildHeight = firstChild.getHeight();
        double secondChildHeight = secondChild.getHeight();

        double minBranchLength = currentPartition.sampleMinBranchLength(vertex);
        double vertexHeight = Math.max(
                firstChildHeight + minBranchLength,
                secondChildHeight + minBranchLength
        );
        vertex.setHeight(vertexHeight);

        return currentClade;
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
        // when a partition has no registered occurrences more, we remove it
        if (partition.getNumberOfOccurrences() <= 0) {
            clade.removePartition(partition);
            return true;
        }
        return false;
    }

    /* -- OTHER METHODS -- */

    public List<BCCDCladePartition> getAllPartitions() {
        Set<BCCDCladePartition> partitions = new HashSet<>();

        for (Clade clade : this.getClades()) {
            for (CladePartition partition : clade.partitions) {
                partitions.add((BCCDCladePartition) partition);
            }
        }

        return new ArrayList<>(partitions);
    }

    @Override
    public String toString() {
        return "CCD1 " + super.toString();
    }

    @Override
    public AbstractCCD copy() {
        CCD1 copy = new CCD1(this.getSizeOfLeavesArray(), false);
        copy.baseTrees.add(this.getSomeBaseTree());
        copy.numBaseTrees = this.getNumberOfBaseTrees();

        AbstractCCD.buildCopy(this, copy);

        return copy;
    }

    @Override
    protected double getNumberOfParameters() {
        return 3 * this.getNumberOfCladePartitions();
    }

}
