package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class BCCDMuSigmaBetaOldOld extends BCCDParameterEstimator {
    @Override
    public void estimateParameters(List<BCCDCladePartition> partitions) {
        double[] solution = new double[2 * partitions.size() + 1];

        updateBeta(partitions, solution);
        updateMus(partitions, solution);
        updateSigmas(partitions, solution);

        double beta = solution[solution.length - 1];

        for (int i = 0; i < partitions.size(); i++) {
            double mu = solution[i];
            double sigma = solution[partitions.size() + i];

            BCCDCladePartition partition = partitions.get(i);
            partition.setLogMeanFunc(x -> mu + beta * x.logBranchLengthOldOld());
            partition.setLogVarianceFunc(x -> sigma);
        }
    }

    public void updateBeta(List<BCCDCladePartition> partitions, double[] parameters) {
        double beta = 0.0;
        int numObservations = 0;

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double[] bs = partition.getObservedLogBranchLengthsOld().toArray();
            double[] bos = partition.getObservedLogBranchLengthsOldOld().toArray();

            if (bs.length < 2) {
                continue;
            }

            double partitionBeta = new Covariance().covariance(bs, bos) / new Variance().evaluate(bos);

            if (Double.isNaN(partitionBeta)) {
                continue;
            }

            beta += bos.length * partitionBeta;
            numObservations += bos.length;
        }

        beta /= numObservations;
        parameters[parameters.length - 1] = beta;
    }

    public void updateMus(List<BCCDCladePartition> partitions, double[] parameters) {
        double beta = parameters[parameters.length - 1];

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double b = partition.getObservedLogBranchLengthsOld().average().orElseThrow();
            double bo = partition.getObservedLogBranchLengthsOldOld().average().orElseThrow();
            parameters[i] = b - beta * bo;
        }
    }

    public void updateSigmas(List<BCCDCladePartition> partitions, double[] parameters) {
        double beta = parameters[parameters.length - 1];

        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Integer> idxWithNegativeEstimate = new ArrayList<>();

        List<Double> allSigmas = new LinkedList<>();

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double[] bs = partition.getObservedLogBranchLengthsOld().toArray();
            double[] bos = partition.getObservedLogBranchLengthsOldOld().toArray();

            double sigma = new Variance().evaluate(bs) - Math.pow(beta, 2) * new Variance().evaluate(bos);

            if (bos.length < 2) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            if (sigma <= 0.0) {
                idxWithNegativeEstimate.add(i);
                continue;
            }

            parameters[partitions.size() + i] = sigma;
            allSigmas.add(sigma);
        }

        double meanSigma = allSigmas.stream().reduce((a, b) -> a + b).get() / partitions.size();
        for (int i : idxWithoutEnoughData) {
            parameters[partitions.size() + i] = meanSigma;
        }

        double lowSigma = new Percentile().evaluate(allSigmas.stream().mapToDouble(x -> x).toArray(), 5);
        for (int i : idxWithNegativeEstimate) {
            parameters[partitions.size() + i] = lowSigma;
        }
    }
}
