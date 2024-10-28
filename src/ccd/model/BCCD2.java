package ccd.model;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer;

import java.util.*;

public class BCCD2 extends BCCD {

    /* -- CONSTRUCTORS & CONSTRUCTION METHODS -- */

    public BCCD2(List<Tree> trees, double burnin) {
        super(trees, burnin);
    }

    public BCCD2(TreeSet treeSet) {
        this(treeSet, false);
    }

    public BCCD2(TreeSet treeSet, boolean storeBaseTrees) {
        super(treeSet, storeBaseTrees);
    }

    public BCCD2(int numLeaves, boolean storeBaseTrees) {
        super(numLeaves, storeBaseTrees);
    }

    @Override
    protected double getNumberOfParameters() {
        return 3 * this.getNumberOfCladePartitions() + 1;
    }

    @Override
    public void initialize() {
        this.estimateParameters();
    }

    @Override
    protected void initializeRootClade(int numLeaves) {
        this.leafArraySize = numLeaves;

        BitSet rootBitSet = BitSet.newBitSet(leafArraySize);
        rootBitSet.set(0, numLeaves);

        this.rootClade = new BCCD2Clade(rootBitSet, this);
        cladeMapping.put(rootClade.getCladeInBits(), rootClade);
    }

    /* -- CLADE CREATION - CLADE CREATION -- */

    @Override
    protected Clade addNewClade(BitSet cladeInBits) {
        Clade clade = new BCCD2Clade(cladeInBits, this);
        cladeMapping.put(cladeInBits, clade);
        return clade;
    }

    /* -- CLADE CREATION - CLADE CREATION -- */

    protected void estimateParameters() {
        List<BCCD2CladePartition> partitions = this.getAllPartitions();

        BCCD2MLE mleProblem = new BCCD2MLE(partitions);
        NonLinearConjugateGradientOptimizer optimizer = new NonLinearConjugateGradientOptimizer(
                NonLinearConjugateGradientOptimizer.Formula.FLETCHER_REEVES,
                new SimpleValueChecker(1e-4, 0)
        );

        PointValuePair solution = optimizer.optimize(
                GoalType.MAXIMIZE,
                new InitialGuess(mleProblem.getInitialGuess()),
                mleProblem.logMLE(),
                mleProblem.logMLEGradient(),
                new MaxEval(50000)
        );

        for (int i = 0; i < partitions.size(); i++) {
            BCCD2CladePartition partition = partitions.get(i);
            partition.setMu(solution.getPoint()[i]);
            partition.setSigma(solution.getPoint()[partitions.size() + i]);
            partition.setBeta(solution.getPoint()[2*partitions.size()]);
        }
    }

    /* -- GENERAL - GENERAL -- */

    public List<BCCD2CladePartition> getAllPartitions() {
        Set<BCCD2CladePartition> partitions = new HashSet<>();

        for (Clade clade : this.getClades()) {
            for (CladePartition partition : clade.partitions) {
                partitions.add((BCCD2CladePartition) partition);
            }
        }

        return new ArrayList<>(partitions);
    }
}
