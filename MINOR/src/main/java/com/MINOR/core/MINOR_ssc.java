package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

/**
 * ssc: simple speed constraint, only check speed constraint after all the iterations.
 */
@Deprecated
public class MINOR_ssc extends MINOR_base {
    protected double f;
    protected double percentile;

    public MINOR_ssc(TimeSeries ts, int p, double threshold, long iterationLimit) {
        this.modelName = "MONOR-ssc";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;
        this.iterationLimit = iterationLimit;
        this.percentile = 99;
        this.f = 1;

        this.isConverge = false;
        this.dim = ts.getDimension();
        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, dim * p);
        V = new DMatrixRMaj(tsLength - p, dim);
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        this.phi = new DMatrixRMaj(p * dim, dim);
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocess();
        MatrixPruning();
        int k = 1;
        this.isConverge = false;
        while (true) {
            estimateIC();
            candidate();
            evaluate();
            if (isConverge || k >= iterationLimit) {
                speedCheck();
                System.out.println("[INFO]Done " + this.modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        double rmse = MyMathUtils.calculateRMSE(ts);
        return new RepairResult(rmse, duration, k);
    }

    protected double[] calculateSpeeds(double f) {
        if (ts.s == null || tsLength < 2 || f <= 0) {
            throw new IllegalArgumentException("[WARNING]Invalid input data when calculate speed.");
        }
        double[] speeds = new double[tsLength];
        int step = (int) Math.max(1, 1 / f);
        for (int i = step; i < tsLength; i += step) {
            DMatrixRMaj p1 = ts.getObsAt(i);
            DMatrixRMaj p2 = ts.getObsAt(i - 1);
            double speed = MyMathUtils.calL2Norm(p2, p1, dim) / Math.abs((ts.getTimestampAt(i) - ts.getTimestampAt(i - 1)));
            speeds[i] = speed;
        }
        return speeds;
    }

    protected double calculatePercentile(double[] speeds, double percentile) {
        Percentile p = new Percentile();
        p.setData(speeds);
        return p.evaluate(percentile);
    }

    /**
     * 对迭代修复后的数据进行速度检查，违约的数据根据情况选择保留或者替换为原始值
     */
    protected void speedCheck() {
        double speedConstraint = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        System.out.println("[INFO]Speed Constraint=" + speedConstraint);
        DMatrixRMaj residual = new DMatrixRMaj(dim, 1);
        DMatrixRMaj result = new DMatrixRMaj(dim, 1);
        for (int t = p; t < tsLength; t++) {
            if (!ts.getLabelAt(t)) {
                double v = MyMathUtils.calL2Norm(ts.getRepairAt(t), ts.getRepairAt(t - 1), dim) / (ts.getTimestampAt(t) - ts.getTimestampAt(t - 1));
                if (v > speedConstraint) {
                    double v_dirty = MyMathUtils.calL2Norm(ts.getObsAt(t), ts.getRepairAt(t - 1), dim) / (ts.getTimestampAt(t) - ts.getTimestampAt(t - 1));
                    if (v > v_dirty) {
                        // 超过约束的值对修复进行缩放,但如果脏值速度更大那就保留修复
//                    double ratio = speedConstraint / v;
                        double ratio = 0;
                        CommonOps_DDRM.subtract(ts.getRepairAt(t), ts.getObsAt(t), residual);
                        CommonOps_DDRM.scale(ratio, residual);
                        CommonOps_DDRM.add(residual, ts.getObsAt(t), result);
                        ts.setRepairAt(t, new DMatrixRMaj(result));
                    }
                }
            }
        }
    }
}
