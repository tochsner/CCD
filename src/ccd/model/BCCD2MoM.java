package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * This class implements the log MLE for the BCCD2 model (up to some constant).
 * <p>
 * Format of the parameter vector: [mu_s ... sigma_s ... beta]
 */
public class BCCD2MoM extends BCCD2MLE {

    public BCCD2MoM(List<BCCD2CladePartition> partitions) {
        super(partitions);
    }

    public void updateBeta(double[] parameters) {
        double beta = 0.0;
        int numObservations = 0;

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);

            double[] bs = partition.getObservedMinLogBranchLength().toArray();
            double[] bds = partition.getObservedMinLogBranchLengthsDown().toArray();

            if (bs.length < 2) {
                continue;
            };

            double partitionBeta = new Covariance().covariance(bs, bds) / new Variance().evaluate(bds);

            if (Double.isNaN(partitionBeta)) {
                continue;
            };

            beta += partitionBeta;
            numObservations += bds.length;
        }

        beta /= numObservations;
        parameters[parameters.length - 1] = beta;
    }

    public void updateMus(double[] parameters) {
        double beta = parameters[parameters.length - 1];

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);

            double b = partition.getObservedMinLogBranchLength().average().orElseThrow();
            double bd  = partition.getObservedMinLogBranchLengthsDown().average().orElseThrow();
            parameters[i] = b - beta*bd;
        }
    }

    public void updateSigmas(double[] parameters) {
        double beta = parameters[parameters.length - 1];

        List<Integer> unsetIndices = new LinkedList<>();

        double totalSigma = 0.0;

        for (int i = 0; i < this.partitions.size(); i++) {
            BCCD2CladePartition partition = this.partitions.get(i);

            double[] bs = partition.getObservedMinLogBranchLength().toArray();
            double[] bds = partition.getObservedMinLogBranchLengthsDown().toArray();

            double sigma = new Variance().evaluate(bs) - Math.pow(beta, 2) * new Variance().evaluate(bds);

            if (bds.length < 2) {
                unsetIndices.add(i);
                continue;
            }

            if (sigma <= 0.0) {
                unsetIndices.add(i);
                continue;
            }

            parameters[this.partitions.size() + i] = sigma;
            totalSigma += sigma;
        }

        double meanSigma = totalSigma / this.partitions.size();
        for (int i : unsetIndices) {
            parameters[this.partitions.size() + i] = meanSigma;
        }
    }

    public double[] estimateParameters() {
        double[] solution = new double[2 * this.partitions.size() + 1];

        this.updateBeta(solution);
        this.updateMus(solution);
        this.updateSigmas(solution);

        return solution;
    }
}
