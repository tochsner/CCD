package ccd.model;

import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class BCCDMuSigmaBetaOldOldOldYoung extends BCCDParameterEstimator {
    @Override
    public void estimateParameters(List<BCCDCladePartition> partitions) {
        double betaOld = getBetaOld(partitions);
        double betaYoung = getbetaYoung(partitions);
        double[] mus = getMus(partitions, betaOld, betaYoung);
        double[] sigmas = getSigmas(partitions, betaOld, betaYoung);

        for (int i = 0; i < partitions.size(); i++) {
            double mu = mus[i];
            double sigma = sigmas[i];

            BCCDCladePartition partition = partitions.get(i);
            partition.setLogMeanFunc(x -> mu + betaOld * x.logBranchLengthOldOld() + betaYoung * x.logBranchLengthOldYoung());
            partition.setLogVarianceFunc(x -> sigma);
        }
    }

    public double getBetaOld(List<BCCDCladePartition> partitions) {
        double beta = 0.0;
        int numObservations = 0;

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double[] bs = partition.getObservedLogBranchLengthsOld().toArray();
            double[] bds = partition.getObservedLogBranchLengthsOldOld().toArray();

            if (bs.length < 2) {
                continue;
            }

            double partitionBeta = new Covariance().covariance(bs, bds) / new Variance().evaluate(bds);

            if (Double.isNaN(partitionBeta)) {
                continue;
            }

            beta += bds.length * partitionBeta;
            numObservations += bds.length;
        }

        beta /= numObservations;
        return beta;
    }

    public double getbetaYoung(List<BCCDCladePartition> partitions) {
        double beta = 0.0;
        int numObservations = 0;

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double[] bs = partition.getObservedLogBranchLengthsOld().toArray();
            double[] bds = partition.getObservedLogBranchLengthsOldYoung().toArray();

            if (bs.length < 2) {
                continue;
            }

            double partitionBeta = new Covariance().covariance(bs, bds) / new Variance().evaluate(bds);

            if (Double.isNaN(partitionBeta)) {
                continue;
            }

            beta += bds.length * partitionBeta;
            numObservations += bds.length;
        }

        beta /= numObservations;
        return beta;
    }

    public double[] getMus(List<BCCDCladePartition> partitions, double betaOld, double betaYoung) {
        double[] mus = new double[partitions.size()];

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double b = partition.getObservedLogBranchLengthsOld().average().orElseThrow();
            double bds = partition.getObservedLogBranchLengthsOldOld().average().orElseThrow();
            double bdl = partition.getObservedLogBranchLengthsOldYoung().average().orElseThrow();
            mus[i] = b - betaOld * bds - betaYoung * bdl;
        }

        return mus;
    }

    public double[] getSigmas(List<BCCDCladePartition> partitions, double betaOld, double betaYoung) {
        double[] sigmas = new double[partitions.size()];

        List<Integer> idxWithoutEnoughData = new ArrayList<>();
        List<Integer> idxWithNegativeEstimate = new ArrayList<>();

        List<Double> allSigmas = new LinkedList<>();

        for (int i = 0; i < partitions.size(); i++) {
            BCCDCladePartition partition = partitions.get(i);

            double[] bs = partition.getObservedLogBranchLengthsOld().toArray();
            double[] bds = partition.getObservedLogBranchLengthsOldOld().toArray();
            double[] bdl = partition.getObservedLogBranchLengthsOldYoung().toArray();

            double sigma = new Variance().evaluate(bs) - Math.pow(betaOld, 2) * new Variance().evaluate(bds) - Math.pow(betaYoung, 2) * new Variance().evaluate(bdl);

            if (bds.length < 2) {
                idxWithoutEnoughData.add(i);
                continue;
            }

            if (sigma <= 0.0) {
                idxWithNegativeEstimate.add(i);
                continue;
            }

            sigmas[i] = sigma;
            allSigmas.add(sigma);
        }

        double meanSigma = allSigmas.stream().reduce((a, b) -> a + b).get() / partitions.size();
        for (int i : idxWithoutEnoughData) {
            sigmas[i] = meanSigma;
        }

        double lowSigma = new Percentile().evaluate(allSigmas.stream().mapToDouble(x -> x).toArray(), 5);
        for (int i : idxWithNegativeEstimate) {
            sigmas[i] = lowSigma;
        }

        return sigmas;
    }
}
