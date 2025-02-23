package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimePoint;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

import java.util.ArrayList;
import java.util.List;

public class MINOR_B_no_iter extends MINOR_U_no_iter {
    protected double[] stRevList;
    protected double[] voRevList;
    protected int firstLabelRev;

    public MINOR_B_no_iter() {

    }

    public MINOR_B_no_iter(TimeSeries ts, int p, double threshold) {
        this.modelName = "MINOR-B-w/o-iter";

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
        this.stRevList = new double[tsLength];
        this.voRevList = new double[tsLength];
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocess();
        reversePreProcess();
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
     * memorize vo^r and st^r
     */
    protected void reversePreProcess() {
        List<TimePoint> rs = new ArrayList<>();
        TimePoint tp;
        int lastLabel = -1;
        int iRev;
        for (int i = tsLength - 1; i >= 0; i--) {
            iRev = tsLength - 1 - i;
            tp = new TimePoint(ts.s.get(i));
            rs.add(tp);
            if (i != tsLength - 1) {
                if (lastLabel == -1) {
                    voRevList[iRev] = MyMathUtils.calL2Norm(tp.getValObs(), rs.get(0).getValRepaired(), dim)
                            / Math.abs(tp.getTimestamp() - rs.get(0).getTimestamp());
                } else {
                    voRevList[iRev] = MyMathUtils.calL2Norm(tp.getValObs(), rs.get(lastLabel).getValTruth(), dim)
                            / Math.abs(tp.getTimestamp() - rs.get(lastLabel).getTimestamp());
                }
            }
            if (tp.getLabel()) {
                if (lastLabel == -1) {
                    firstLabelRev = iRev;
                    stRevList[iRev] = MyMathUtils.calL2Norm(tp.getValTruth(), rs.get(0).getValRepaired(), dim)
                            / Math.abs(tp.getTimestamp() - rs.get(0).getTimestamp());
                } else {
                    stRevList[iRev] = MyMathUtils.calL2Norm(tp.getValTruth(), rs.get(lastLabel).getValTruth(), dim)
                            / Math.abs(tp.getTimestamp() - rs.get(lastLabel).getTimestamp());
                }
                lastLabel = iRev;
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
    @Override
    protected double validate(int t) {
        double s_t = Double.MAX_VALUE;
        double s_t_r = Double.MAX_VALUE;
        double v;
        double vRev;
        int lastLabelP = -1;
        int nextLabelP = -1;
        TimeSeries tsUsed = this.ts;
        if (t == 0 || t == tsLength - 1) {
            // the first and last data points will not be modified.
            return 0;
        }
        // Find the two preceding markers and calculate the speed to use as the speed constraint s_t
        for (int i = t - 1; i >= 0; i--) {
            if (i == 0) {
                lastLabelP = 0;
                s_t = this.s_g;
                break;
            }
            if (tsUsed.getLabelAt(i)) {
                // Find the nearest preceding markers.
                lastLabelP = i;
                s_t = stList[i];
                break;
            }
        }
        // Find the two subsequent markers and calculate the speed to use as the speed constraint s_t_r
        for (int i = t + 1; i < tsLength; i++) {
            if (i == tsLength - 1) {
                nextLabelP = tsLength - 1;
                s_t_r = s_g;
                break;
            }
            if (tsUsed.getLabelAt(i)) {
                nextLabelP = i;
                s_t_r = stRevList[tsLength - 1 - i];
                break;
            }
        }
        vRev = MyMathUtils.calL2Norm(tsUsed.getCandidateAt(t), tsUsed.getRepairAt(nextLabelP), dim)
                / Math.abs((tsUsed.getTimestampAt(t) - tsUsed.getTimestampAt(nextLabelP)));

        v = MyMathUtils.calL2Norm(tsUsed.getCandidateAt(t), tsUsed.getRepairAt(lastLabelP), dim)
                / Math.abs((tsUsed.getTimestampAt(t) - tsUsed.getTimestampAt(lastLabelP)));

            // compute validity from both directions
            double minSC = Math.min(Math.min(s_g, s_t), s_t_r);
            if (v > minSC || vRev > minSC) {
                double vo, voRev;
                double validity;
                vo = voList[t];
                voRev = voRevList[tsLength - 1 - t];
                double validityRev;
                double invalidity1 = v - s_t;
                double invalidity2 = vo - s_t;
                if (invalidity1 < 0) {
                    // The candidate satisfies the speed constraint on this side.
                    validity = 1;
                } else {
                    if (invalidity2 < 0) {
                        // The candidate does not satisfy the speed constraint on this side,
                        // but the observe value met the speed constraint, so validity of this side is negative.
                        validity = Math.max(-1, 1 + invalidity1 / invalidity2);
                    } else {
                        validity = Math.min(1, 1 - invalidity1 / invalidity2);
                    }
                }
                double invalidity1Rev = vRev - s_t_r;
                double invalidity2Rev = voRev - s_t_r;
                if (invalidity1Rev < 0) {
                    validityRev = 1;
                } else {
                    if (invalidity2Rev < 0) {
                        validityRev = Math.max(-1, 1 + invalidity1Rev / invalidity2Rev);
                    } else {
                        validityRev = Math.min(1, 1 - invalidity1Rev / invalidity2Rev);
                    }
                }
                long dist = Math.abs(tsUsed.getTimestampAt(t) - tsUsed.getTimestampAt(lastLabelP));
                long distRev = Math.abs(tsUsed.getTimestampAt(nextLabelP) - tsUsed.getTimestampAt(t));
                return (validity * Math.exp(distRev) + validityRev * Math.exp(dist)) / (Math.exp(dist) + Math.exp(distRev));
            }
            return 1;
    }
}
