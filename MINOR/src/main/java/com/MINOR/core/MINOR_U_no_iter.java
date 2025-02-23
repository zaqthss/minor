package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

public class MINOR_U_no_iter extends VARX {
    protected double s_g;   // global speed constraint
    protected double percentile;
    protected double f;
    protected int fistLabel;
    protected double[] stList;
    protected double[] voList;


    public MINOR_U_no_iter() {

    }

    public MINOR_U_no_iter(TimeSeries ts, int p, double threshold) {
        this.modelName = "MINOR-U-w/o-iter";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;

        this.percentile = 99;
        this.f = 1;

        this.dim = ts.getDimension();
        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, dim * p);
        V = new DMatrixRMaj(tsLength - p, dim);
        this.phi = new DMatrixRMaj(p * dim, dim);

        this.s_g = 0;
        this.stList = new double[tsLength];
        this.voList = new double[tsLength];
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocess();
        OLSPre();
        OLS();
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        repair();
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]" + this.modelName + "(" + p + ") Done.");
        return new RepairResult(MyMathUtils.calculateRMSE(ts), duration, 1);
    }

    /**
     * initialize data and memorize vo and st
     */
    @Override
    protected void preprocess() {
        int lastLabel = -1;
        boolean isLabeled;
        for (int i = 0; i < tsLength; i++) {
            isLabeled = false;
            DMatrixRMaj mat;
            DMatrixRMaj obs = ts.getObsAt(i);
            DMatrixRMaj truth = ts.getTruthAt(i);
            if (ts.getLabelAt(i)) {
                mat = new DMatrixRMaj(truth);
                isLabeled = true;
            } else {
                mat = new DMatrixRMaj(obs);
            }
            ts.setRepairAt(i, new DMatrixRMaj(mat));
            ts.setLastRepairAt(i, new DMatrixRMaj(mat));
            ts.setLastModifyAt(i, new DMatrixRMaj(mat));
            if (isLabeled) {
                if (lastLabel != -1) {
                    // The speed between labeled data points will not change after preprocessing.
                    stList[i] = MyMathUtils.calL2Norm(truth, ts.getTruthAt(lastLabel), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(lastLabel));
                } else {
                    this.fistLabel = i;
                    stList[i] = MyMathUtils.calL2Norm(truth, ts.getRepairAt(0), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(0));
                }
                lastLabel = i;
            }
            if (i > 0) {
                if (lastLabel == -1) {
                    voList[i] = MyMathUtils.calL2Norm(obs, ts.getRepairAt(0), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(0));
                } else {
                    voList[i] = MyMathUtils.calL2Norm(obs, ts.getTruthAt(lastLabel), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(lastLabel));
                }
            }
        }
    }

    protected void repair() {
        DMatrixRMaj temp = new DMatrixRMaj(dim, 1);
        DMatrixRMaj phiSubMat = new DMatrixRMaj(dim, dim);
        DMatrixRMaj z = new DMatrixRMaj(dim, 1);
        DMatrixRMaj z2 = new DMatrixRMaj(dim, 1);
        for (int t = p; t < tsLength; t++) {
            if (!ts.getLabelAt(t)) {
                temp.zero();
                for (int i = 1; i <= p; i++) {
                    // Obtain the corresponding submatrix when the lag coefficient is p
                    CommonOps_DDRM.extract(phi, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, phiSubMat);
                    CommonOps_DDRM.subtract(ts.getRepairAt(t - i), ts.getObsAt(t - i), z);
                    CommonOps_DDRM.multTransA(phiSubMat, z, z2);
                    CommonOps_DDRM.add(temp, z2, temp);
                }
                CommonOps_DDRM.add(temp, ts.getObsAt(t), temp);
                // The L2 norm between the predicted value and the original predicted value is greater than the threshold
                double l2norm = MyMathUtils.calL2Norm(temp, ts.getRepairAt(t), dim);
                if (l2norm > this.threshold) {
                    ts.setCandidateAt(t, new DMatrixRMaj(temp));
                    if (validate(t) > 0) {
                        ts.setRepairAt(t, new DMatrixRMaj(temp));
                    }
                }
            }
        }
        ts.clearCandidates();
    }

    /**
     * Check whether a candidate meets the speed constraint s_t and s_g, and provide its validity
     *
     * @param t timesteamp
     * @return validity
     */
    protected double validate(int t) {
        double s_t = Double.MAX_VALUE;
        double v;
        int lastLabelP = 0;
        for (int i = t - 1; i >= 0; i--) {
            if (i == 0) {
                s_t = s_g;
                break;
            }
            if (ts.getLabelAt(i)) {
                // Find the nearest preceding labeled point
                lastLabelP = i;
                s_t = stList[i];
                break;
            }
        }
        v = MyMathUtils.calL2Norm(ts.getCandidateAt(t), ts.getRepairAt(lastLabelP), dim) / (ts.getTimestampAt(t) - ts.getTimestampAt(lastLabelP));
        if (v > s_g || v > s_t) {
            double vo = voList[t];
            if (vo <= s_t) {
                return 0;
            } else {
                double invalidity1 = v - s_t;
                double invalidity2 = vo - s_t;
                return Math.max(1 - invalidity1 / invalidity2, 0);  // validity
            }
        }
        return 1;
    }

    /**
     * Sample and calculate the speed from observed data
     *
     * @param f Sampling ratio
     * @return Speed of the sampled data
     */
    protected double[] calculateSpeeds(double f) {
        if (ts.s == null || tsLength < 2 || f <= 0) {
            throw new IllegalArgumentException("[WARNING]Invalid input data when calculate speed.");
        }
        double[] speeds = new double[tsLength];
        int step = (int) Math.max(1, 1 / f);  // Calculate the step size for each sampling based on the sampling ratio f
        for (int i = step; i < tsLength; i += step) {
            // Calculate the speed of each point relative to the previous adjacent point.
            DMatrixRMaj p1 = ts.getObsAt(i);
            DMatrixRMaj p2 = ts.getObsAt(i - 1);
            // Use the L2 norm as the speed between adjacent points.
            double speed = MyMathUtils.calL2Norm(p2, p1, dim) / Math.abs((ts.getTimestampAt(i) - ts.getTimestampAt(i - 1)));
            speeds[i] = speed;
        }
        return speeds;
    }

    /**
     * Calculate the speed at a given percentile to use as a speed constraint.
     *
     * @param speeds     A sequence of speed values
     * @param percentile The percentile, where 0 < percentile <= 100; 99 is commonly used
     * @return The speed value corresponding to the given percentile
     */
    protected double calculatePercentile(double[] speeds, double percentile) {
        Percentile p = new Percentile();
        p.setData(speeds);
        return p.evaluate(percentile);
    }
}
