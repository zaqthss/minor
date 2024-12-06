package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimePoint;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.ops.SimpleOperations_DDRM;

import java.util.*;

public class MINOR_B extends MINOR_U {
    protected DMatrixRMaj phiReverse;
    protected TimeSeries tsReverse;
    protected DMatrixRMaj ZRev;
    protected DMatrixRMaj VRev;
    protected DMatrixRMaj ARev;
    protected DMatrixRMaj BRev;
    protected DMatrixRMaj XRev;
    protected DMatrixRMaj[] XRevList;

    protected double[] stRevList;
    protected double[] voRevList;

    protected int firstLabelRev;

    protected List<Double> l2normList;

    /**
     * 记录在正向修复和反向修复中同时存在的candidates
     */
    protected Set<Integer> RepairedPointSet;

    public MINOR_B() {

    }


    public MINOR_B(TimeSeries ts, int p, double threshold, long iterationLimit) {
        this.modelName = "MINOR-B";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;
        this.iterationLimit = iterationLimit;
        this.percentile = 99;
        this.f = 1;

        r = -1;

        this.dim = ts.getDimension();
        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, dim * p);
        V = new DMatrixRMaj(tsLength - p, dim);
        X = new DMatrixRMaj(dim, tsLength - p);
        XList = new DMatrixRMaj[tsLength - p];
        ZRev = new DMatrixRMaj(tsLength - p, dim * p);
        VRev = new DMatrixRMaj(tsLength - p, dim);
        XRev = new DMatrixRMaj(dim, tsLength - p);
        XRevList = new DMatrixRMaj[tsLength - p];
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        ARev = new DMatrixRMaj();
        BRev = new DMatrixRMaj();

        this.phi = new DMatrixRMaj(p * dim, dim);
        this.phiReverse = new DMatrixRMaj(p * dim, dim);

        this.RepairedPointSet = new HashSet<>();

        this.s_g = 0;
        this.stList = new double[tsLength];
        this.voList = new double[tsLength];
        this.stRevList = new double[tsLength];
        this.voRevList = new double[tsLength];

        this.l2normList = new ArrayList<>();
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocessIC();
        tsReverse = reversePreProcessIC();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        for (int i = 0; i < tsLength; i++) {
            if (this.ts.getLabelAt(i)) {
                RepairedPointSet.add(i);
            }
        }
        while (true) {
            estimateIC();
            estimateICReverse();
            candidateIC(this.ts, this.phi, this.XList, this.Z);
            candidateIC(this.tsReverse, this.phiReverse, this.XRevList, this.ZRev);
            evaluate();
            updateZ();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + this.modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            if (k % (iterationLimit / 10) == 0) {
                System.out.println(k);
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]runtime: " + duration + " ms");
        double rmse = MyMathUtils.calculateRMSE(ts);
        return new RepairResult(rmse, duration, k);
    }

    /**
     * when running application experiment
     *
     * @return repaired time series
     */
    public TimeSeries runApp() {
        preprocessIC();
        tsReverse = reversePreProcessIC();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        for (int i = 0; i < tsLength; i++) {
            if (this.ts.getLabelAt(i)) {
                RepairedPointSet.add(i);
            }
        }
        while (true) {
            estimateIC();
            estimateICReverse();
            candidateIC(this.ts, this.phi, this.XList, this.Z);
            candidateIC(this.tsReverse, this.phiReverse, this.XList, this.ZRev);
            evaluate();
            updateZ();
            if (isConverge || k >= iterationLimit) {
                break;
            }
            k++;
        }
        return ts;
    }

    @Override
    public RepairResult runWithTimer() {
        long startTime = System.nanoTime();
        preprocessIC();
        tsReverse = reversePreProcessIC();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        double[] timeCost = new double[4];
        for (int i = 0; i < tsLength; i++) {
            if (this.ts.getLabelAt(i)) {
                RepairedPointSet.add(i);
            }
        }
        while (true) {
            long startTime1 = System.nanoTime();
            estimateIC();
            estimateICReverse();
            long endTime1 = System.nanoTime();
            candidateIC(this.ts, this.phi, this.XList, this.Z);
            candidateIC(this.tsReverse, this.phiReverse, this.XRevList, this.ZRev);
            long endTime2 = System.nanoTime();
            evaluate();
            long endTime3 = System.nanoTime();
            updateZ();
            long endTime4 = System.nanoTime();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + this.modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            if (k % (iterationLimit / 10) == 0) {
                System.out.println("[REPAIRING] " + k + " / " + iterationLimit);
            }
            k++;
            timeCost[0] += ((double) (endTime1 - startTime1)) / 1_000_000.0d;
            timeCost[1] += ((double) (endTime2 - endTime1)) / 1_000_000.0d;
            timeCost[2] += ((double) (endTime3 - endTime2)) / 1_000_000.0d;
            timeCost[3] += ((double) (endTime4 - endTime3)) / 1_000_000.0d;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]runtime: " + duration + " ms");
        double rmse = MyMathUtils.calculateRMSE(ts);
        RepairResult result = new RepairResult(rmse, duration, k);
        System.out.println(result);
        System.out.println("[INFO]time cost of 4 parts: estimate: " + timeCost[0]
                + ", candidate: " + timeCost[1]
                + ", evaluate: " + timeCost[2]
                + ", updateZ: " + timeCost[3]);
        return result;
    }

    @Override
    public RepairResult runSD() {
        long startTime = System.nanoTime();
        preprocess();
        tsReverse = reversePreProcess();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        double[] timeCost = new double[3];
        for (int i = 0; i < tsLength; i++) {
            if (this.ts.getLabelAt(i)) {
                RepairedPointSet.add(i);
            }
        }
        while (true) {
            long startTime1 = System.nanoTime();
            estimate();
            estimateReverse();
            long endTime1 = System.nanoTime();
            candidate(this.ts, this.phi);
            candidate(this.tsReverse, this.phiReverse);
            long endTime2 = System.nanoTime();
            evaluate();
            long endTime3 = System.nanoTime();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + this.modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            if (k % (iterationLimit / 10) == 0) {
                System.out.println("[REPAIRING] " + k + " / " + iterationLimit);
            }
            k++;
            timeCost[0] += ((double) (endTime1 - startTime1)) / 1_000_000.0d;
            timeCost[1] += ((double) (endTime2 - endTime1)) / 1_000_000.0d;
            timeCost[2] += ((double) (endTime3 - endTime2)) / 1_000_000.0d;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]runtime: " + duration + " ms");
        double rmse = MyMathUtils.calculateRMSE(ts);
        RepairResult result = new RepairResult(rmse, duration, k);
        System.out.println(result);
        System.out.println("[INFO]time cost of 3 parts: estimate: " + timeCost[0]
                + ", candidate: " + timeCost[1]
                + ", evaluate: " + timeCost[2]);
        return result;
    }

    @Override
    public RepairResult runWithoutCandiIC() {
        long startTime = System.nanoTime();
        preprocessIC();
        tsReverse = reversePreProcess();
        MatrixPruning();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        double[] timeCost = new double[3];
        for (int i = 0; i < tsLength; i++) {
            if (this.ts.getLabelAt(i)) {
                RepairedPointSet.add(i);
            }
        }
        while (true) {
            long startTime1 = System.nanoTime();
            estimateIC();
            estimateICReverse();
            long endTime1 = System.nanoTime();
            candidate(this.ts, this.phi);
            candidate(this.tsReverse, this.phiReverse);
            long endTime2 = System.nanoTime();
            evaluate();
            long endTime3 = System.nanoTime();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + this.modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            if (k % (iterationLimit / 10) == 0) {
                System.out.println("[REPAIRING] " + k + " / " + iterationLimit);
            }
            k++;
            timeCost[0] += ((double) (endTime1 - startTime1)) / 1_000_000.0d;
            timeCost[1] += ((double) (endTime2 - endTime1)) / 1_000_000.0d;
            timeCost[2] += ((double) (endTime3 - endTime2)) / 1_000_000.0d;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]runtime: " + duration + " ms");
        double rmse = MyMathUtils.calculateRMSE(ts);
        RepairResult result = new RepairResult(rmse, duration, k);
        System.out.println(result);
        System.out.println("[INFO]time cost of 3 parts: estimate: " + timeCost[0]
                + ", candidate: " + timeCost[1]
                + ", evaluate: " + timeCost[2]);
        return result;
    }

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

    @Override
    protected void preprocessIC() {
        int lastLabel = -1;
        boolean isLabeled;
        List<DMatrixRMaj> z = new ArrayList<>();
        X = new DMatrixRMaj(dim, tsLength - p);
        for (int i = 0; i < tsLength; i++) {
            isLabeled = false;
            DMatrixRMaj mat;
            DMatrixRMaj obs = ts.getObsAt(i);
            DMatrixRMaj truth = ts.getTruthAt(i);
            if (i >= p) {
                for (int j = 0; j < dim; j++) {
                    // init obs matrix X
                    X.set(j, i - p, obs.get(j, 0));
                }
            }
            if (ts.getLabelAt(i)) {
                mat = new DMatrixRMaj(truth);
                isLabeled = true;
                z.add(CommonOps_DDRM.subtract(truth, obs, null));
            } else {
                mat = new DMatrixRMaj(obs);
                z.add(new DMatrixRMaj(dim, 1));
            }
            ts.setRepairAt(i, new DMatrixRMaj(mat));
            ts.setLastRepairAt(i, new DMatrixRMaj(mat));
            ts.setLastModifyAt(i, new DMatrixRMaj(mat));
            if (isLabeled) {
                if (lastLabel != -1) {
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
        CommonOps_DDRM.columnsToVector(X, XList);
        // matrix Z
        for (int i = 0; i < tsLength - p; i++) {
            for (int j = 0; j < p; j++) {
                for (int k = 0; k < dim; k++) {
                    Z.set(i, dim * j + k, z.get(p - 1 + i - j).get(k, 0));
                }
            }
        }
        for (int i = 0; i < tsLength - p; i++) {
            for (int k = 0; k < dim; k++) {
                V.set(i, k, z.get(p + i).get(k, 0));
            }
        }
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        DMatrixRMaj[] wrapper = {Z, V};
        MyMathUtils.MatrixPruning(wrapper);
        A.reshape(Z.numCols, Z.numCols);
        B.reshape(Z.numCols, V.numCols);
        CommonOps_DDRM.multTransA(wrapper[0], wrapper[0], A);         // A_0 = Z'*Z
        CommonOps_DDRM.multTransA(wrapper[0], wrapper[1], B);
    }

    /**
     * initialize data of the reversed time series and memorize vo^r and st^r
     *
     * @return reversed time series
     */
    protected TimeSeries reversePreProcess() {
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
        return new TimeSeries(rs);
    }

    /**
     * initialize data of the reversed time series ,memorize vo^r and st^r
     * and initialize variables for Incremental Candidate Generation
     *
     * @return reversed time series
     */
    protected TimeSeries reversePreProcessIC() {
        List<TimePoint> rs = new ArrayList<>();
        TimePoint tp;
        List<DMatrixRMaj> zR = new ArrayList<>();
        int lastLabel = -1;
        int iRev;
        for (int i = tsLength - 1; i >= 0; i--) {
            iRev = tsLength - 1 - i;
            tp = new TimePoint(ts.s.get(i));
            rs.add(tp);
            if (iRev >= p) {
                for (int j = 0; j < dim; j++) {
                    // init obs matrix X
                    XRev.set(j, iRev - p, tp.getValObs().get(j, 0));
                }
            }
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
                zR.add(CommonOps_DDRM.subtract(tp.getValTruth(), tp.getValObs(), null));
                lastLabel = iRev;
            } else {
                zR.add(new DMatrixRMaj(dim, 1));
            }
        }
        CommonOps_DDRM.columnsToVector(XRev, XRevList);
        for (int i = 0; i < tsLength - p; i++) {
            for (int j = 0; j < p; j++) {
                for (int k = 0; k < dim; k++) {
                    ZRev.set(i, dim * j + k, zR.get(p - 1 + i - j).get(k, 0));
                }
            }
        }
        for (int i = 0; i < tsLength - p; i++) {
            for (int k = 0; k < dim; k++) {
                VRev.set(i, k, zR.get(p + i).get(k, 0));
            }
        }
        DMatrixRMaj[] wrapperReverse = {ZRev, VRev};
        MyMathUtils.MatrixPruning(wrapperReverse);
        A.reshape(Z.numCols, Z.numCols);
        B.reshape(Z.numCols, V.numCols);
        ARev.reshape(ZRev.numCols, ZRev.numCols);
        BRev.reshape(ZRev.numCols, VRev.numCols);
        CommonOps_DDRM.multTransA(wrapperReverse[0], wrapperReverse[0], ARev);
        CommonOps_DDRM.multTransA(wrapperReverse[0], wrapperReverse[1], BRev);
        return new TimeSeries(rs);
    }

    @Override
    protected void OLSPre() {
        List<DMatrixRMaj> z = new ArrayList<>();
        List<DMatrixRMaj> zR;
        for (int i = 0; i < tsLength; i++) {
            DMatrixRMaj resultMatrix = new DMatrixRMaj(dim, 1);
            CommonOps_DDRM.subtract(ts.getRepairAt(i), ts.getObsAt(i), resultMatrix);
            z.add(resultMatrix);
        }
        zR = new ArrayList<>(z);
        Collections.reverse(zR);
        // matrix Z
        for (int i = 0; i < tsLength - p; i++) {
            for (int j = 0; j < p; j++) {
                for (int k = 0; k < dim; k++) {
                    Z.set(i, dim * j + k, z.get(p - 1 + i - j).get(k, 0));
                    ZRev.set(i, dim * j + k, zR.get(p - 1 + i - j).get(k, 0));
                }
            }
        }
        // matrix V
        for (int i = 0; i < tsLength - p; i++) {
            for (int k = 0; k < dim; k++) {
                V.set(i, k, z.get(p + i).get(k, 0));
                VRev.set(i, k, zR.get(p + i).get(k, 0));
            }
        }
    }

    @Override
    protected void MatrixPruning() {
        OLSPre();
        DMatrixRMaj[] wrapper = {Z, V};
        DMatrixRMaj[] wrapperReverse = {ZRev, VRev};
        MyMathUtils.MatrixPruning(wrapper);
        MyMathUtils.MatrixPruning(wrapperReverse);
        Z = wrapper[0];
        V = wrapper[1];
        ZRev = wrapperReverse[0];
        VRev = wrapperReverse[1];
        A.reshape(Z.numCols, Z.numCols);
        B.reshape(Z.numCols, V.numCols);
        ARev.reshape(ZRev.numCols, ZRev.numCols);
        BRev.reshape(ZRev.numCols, VRev.numCols);
        CommonOps_DDRM.multTransA(Z, Z, A);         // A_0 = Z'*Z
        CommonOps_DDRM.multTransA(Z, V, B);         // B_0 = Z'*V
        CommonOps_DDRM.multTransA(ZRev, ZRev, ARev);
        CommonOps_DDRM.multTransA(ZRev, VRev, BRev);
    }

    /**
     *  learn parameter Phi^r without Matrix Pruning and Incremental Computation
     */
    protected void estimateReverse() {
        OLSPre();
        OLSReverse();
    }

    /**
     * performing OLS for Phi^r of reversed time series
     */
    protected void OLSReverse() {
        DMatrixRMaj result_temp1 = new DMatrixRMaj(dim * p, dim * p);
        DMatrixRMaj result_temp2 = new DMatrixRMaj(dim * p, tsLength - p);
        DMatrixRMaj result_invert = new DMatrixRMaj(dim * p, dim * p);
        CommonOps_DDRM.multTransA(ZRev, ZRev, result_temp1);
        CommonOps_DDRM.invert(result_temp1, result_invert);
        if (Double.isInfinite(result_invert.get(0, 0)) || Double.isNaN(result_invert.get(0, 0))) {
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(result_temp1, result_invert);
        }
        CommonOps_DDRM.multTransB(result_invert, ZRev, result_temp2);
        CommonOps_DDRM.mult(result_temp2, VRev, phiReverse);
    }

    /**
     * learn parameter Phi^r of reverse repair with Matrix Pruning only
     */
    protected void estimateMRReverse() {
        DMatrixRMaj A_inv = new DMatrixRMaj(dim * p, dim * p);
        // (A)^(-1)
        CommonOps_DDRM.invert(ARev, A_inv);
        if (Double.isInfinite(A_inv.get(0, 0)) || Double.isNaN(A_inv.get(0, 0))) {
            // for singular matrix A, use pseudo-invert instead.
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(ARev, A_inv);
        }
        CommonOps_DDRM.mult(A_inv, BRev, phiReverse);
    }

    /**
     * learn parameter Phi^r of reverse repair with Incremental Computation
     */
    protected void estimateICReverse() {
        if (-1 != r) {
            // 增量式更新矩阵元素
            DMatrixRMaj subMat = new DMatrixRMaj(dim, dim);
            DMatrixRMaj mat1 = new DMatrixRMaj(dim, dim);
            DMatrixRMaj mat2 = new DMatrixRMaj(dim, dim);
            DMatrixRMaj mat3 = new DMatrixRMaj(dim, dim);
            DMatrixRMaj mat4 = new DMatrixRMaj(dim, dim);
            r = tsLength - r - 1;
            r += 1;                   // 让r的取值范围与论文中的r对齐
            for (int i = 1; i <= p; i++) {
                // a_ii
                if (r >= p + 1 - i && r <= tsLength - i) {
                    CommonOps_DDRM.extract(ARev, (i - 1) * dim, (i - 1) * dim + dim, (i - 1) * dim, (i - 1) * dim + dim, subMat);
                    CommonOps_DDRM.multTransB(tsReverse.getResidualAt(r - 1), tsReverse.getResidualAt(r - 1), mat1);
                    CommonOps_DDRM.multTransB(tsReverse.getLastResidualAt(r - 1), tsReverse.getLastResidualAt(r - 1), mat2);
                    CommonOps_DDRM.subtract(mat1, mat2, mat3);
                    CommonOps_DDRM.add(subMat, mat3, mat1);
                    for (int m = 0; m < dim; m++) {
                        for (int n = 0; n < dim; n++) {
                            ARev.set((i - 1) * dim + m, (i - 1) * dim + n, mat1.get(m, n));
                        }
                    }

                }
                // a_ij, a_ji, assert i<j
                for (int j = i + 1; j <= p; j++) {
                    if (r >= p + 1 - j && r < p + 1 - i) {
                        CommonOps_DDRM.extract(ARev, (i - 1) * dim, (i - 1) * dim + dim, (j - 1) * dim, (j - 1) * dim + dim, subMat);
                        CommonOps_DDRM.subtract(tsReverse.getResidualAt(r - 1), tsReverse.getLastResidualAt(r - 1), mat1);
                        CommonOps_DDRM.multTransB(tsReverse.getLastResidualAt(r + j - i - 1), mat1, mat2);
                        CommonOps_DDRM.add(subMat, mat2, mat1);
                        for (int m = 0; m < dim; m++) {
                            for (int n = 0; n < dim; n++) {
                                ARev.set((i - 1) * dim + m, (j - 1) * dim + n, mat1.get(m, n));
                                ARev.set((j - 1) * dim + m, (i - 1) * dim + n, mat1.get(n, m));
                            }
                        }
                    } else if (r > tsLength - j && r <= tsLength - i) {
                        CommonOps_DDRM.extract(ARev, (i - 1) * dim, (i - 1) * dim + dim, (j - 1) * dim, (j - 1) * dim + dim, subMat);
                        CommonOps_DDRM.subtract(tsReverse.getResidualAt(r - 1), tsReverse.getLastResidualAt(r - 1), mat1);
                        CommonOps_DDRM.multTransB(mat1, tsReverse.getLastResidualAt(r - j + i - 1), mat2);
                        CommonOps_DDRM.add(subMat, mat2, mat1);
                        for (int m = 0; m < dim; m++) {
                            for (int n = 0; n < dim; n++) {
                                ARev.set((i - 1) * dim + m, (j - 1) * dim + n, mat1.get(m, n));
                                ARev.set((j - 1) * dim + m, (i - 1) * dim + n, mat1.get(n, m));
                            }
                        }
                    } else if (r >= p + 1 - i && r <= tsLength - j) {
                        CommonOps_DDRM.extract(ARev, (i - 1) * dim, (i - 1) * dim + dim, (j - 1) * dim, (j - 1) * dim + dim, subMat);
                        CommonOps_DDRM.subtract(tsReverse.getResidualAt(r - 1), tsReverse.getLastResidualAt(r - 1), mat1);
                        CommonOps_DDRM.multTransB(tsReverse.getLastResidualAt(r + j - i - 1), mat1, mat2);
                        CommonOps_DDRM.multTransB(mat1, tsReverse.getLastResidualAt(r - j + i - 1), mat3);
                        CommonOps_DDRM.add(mat2, mat3, mat4);
                        CommonOps_DDRM.add(subMat, mat4, mat1);
                        for (int m = 0; m < dim; m++) {
                            for (int n = 0; n < dim; n++) {
                                ARev.set((i - 1) * dim + m, (j - 1) * dim + n, mat1.get(m, n));
                                ARev.set((j - 1) * dim + m, (i - 1) * dim + n, mat1.get(n, m));
                            }
                        }
                    } // else r<p+1-j || r>n-i, no change
                }
                // bi
                if (r >= p + 1 - i && r < p + 1) {
                    CommonOps_DDRM.extract(BRev, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, subMat);
                    CommonOps_DDRM.subtract(tsReverse.getResidualAt(r - 1), tsReverse.getLastResidualAt(r - 1), mat1);
                    CommonOps_DDRM.multTransB(mat1, tsReverse.getLastResidualAt(r + i - 1), mat2);
                    CommonOps_DDRM.add(subMat, mat2, mat1);
                    for (int m = 0; m < dim; m++) {
                        for (int n = 0; n < dim; n++) {
                            BRev.set((i - 1) * dim + m, n, mat1.get(m, n));
                        }
                    }
                } else if (r > tsLength - i) {
                    CommonOps_DDRM.extract(BRev, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, subMat);
                    CommonOps_DDRM.subtract(tsReverse.getResidualAt(r - 1), tsReverse.getLastResidualAt(r - 1), mat1);
                    CommonOps_DDRM.multTransB(tsReverse.getLastResidualAt(r - i - 1), mat1, mat2);
                    CommonOps_DDRM.add(subMat, mat2, mat1);
                    for (int m = 0; m < dim; m++) {
                        for (int n = 0; n < dim; n++) {
                            BRev.set((i - 1) * dim + m, n, mat1.get(m, n));
                        }
                    }
                } else if (r >= p + 1 && r <= tsLength - i) {
                    CommonOps_DDRM.extract(BRev, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, subMat);
                    CommonOps_DDRM.subtract(tsReverse.getResidualAt(r - 1), tsReverse.getLastResidualAt(r - 1), mat1);
                    CommonOps_DDRM.multTransB(mat1, tsReverse.getLastResidualAt(r + i - 1), mat2);
                    CommonOps_DDRM.multTransB(tsReverse.getLastResidualAt(r - i - 1), mat1, mat3);
                    CommonOps_DDRM.add(mat2, mat3, mat4);
                    CommonOps_DDRM.add(subMat, mat4, mat1);
                    for (int m = 0; m < dim; m++) {
                        for (int n = 0; n < dim; n++) {
                            BRev.set((i - 1) * dim + m, n, mat1.get(m, n));
                        }
                    }
                } // else r>p+1-i, no change
            }
            tsReverse.setLastRepairAt(r - 1, new DMatrixRMaj(tsReverse.getRepairAt(r - 1)));
            r -= 1;
            r = tsLength - r - 1;
        }
        DMatrixRMaj A_inv = new DMatrixRMaj(dim * p, dim * p);
        CommonOps_DDRM.invert(ARev, A_inv);
        if (Double.isInfinite(A_inv.get(0, 0)) || Double.isNaN(A_inv.get(0, 0))) {
            // for singular matrix ARev, use pseudo-invert instead.
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(ARev, A_inv);
        }
        CommonOps_DDRM.mult(A_inv, BRev, phiReverse);
    }

    /**
     * generate candidates without Incremental Computation of Z and Z^r
     */
    protected void candidate(TimeSeries TS, DMatrixRMaj Phi) {
        DMatrixRMaj temp = new DMatrixRMaj(dim, 1);
        DMatrixRMaj phiSubMat = new DMatrixRMaj(dim, dim);
        DMatrixRMaj z = new DMatrixRMaj(dim, 1);
        DMatrixRMaj z2 = new DMatrixRMaj(dim, 1);
        for (int t = p; t < tsLength; t++) {
            if (!TS.getLabelAt(t)) {
                temp.zero();
                for (int i = 1; i <= p; i++) {
                    // 获取当滞后系数为p时对应的子矩阵
                    CommonOps_DDRM.extract(Phi, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, phiSubMat);
                    CommonOps_DDRM.subtract(TS.getRepairAt(t - i), TS.getObsAt(t - i), z);
                    CommonOps_DDRM.multTransA(phiSubMat, z, z2);
                    CommonOps_DDRM.add(temp, z2, temp);
                }
                CommonOps_DDRM.add(temp, TS.getObsAt(t), temp);
                // 预测值与原预测值的L2范数大于阈值
                double l2norm = MyMathUtils.calL2Norm(temp, TS.getRepairAt(t), dim);
                if (l2norm > this.threshold) {
                    TS.setCandidateAt(t, new DMatrixRMaj(temp));
                }
            }
        }
    }

    /**
     * generate candidates with Incremental Computation of Z and Z^r
     */
    protected void candidateIC(TimeSeries TS, DMatrixRMaj Phi, DMatrixRMaj[] _XList, DMatrixRMaj _Z) {
        // Y_hat = X + Phi^T * Z^T
        // each column of Y_hat is a candidate at the index
        DMatrixRMaj[] candidates = new DMatrixRMaj[tsLength - p];
        DMatrixRMaj candidate = new DMatrixRMaj(dim, 1);
        CommonOps_DDRM.columnsToVector(CommonOps_DDRM.multTransAB(Phi, _Z, null), candidates);
        for (int t = p; t < tsLength; t++) {
            if (!TS.getLabelAt(t) && !MatrixFeatures_DDRM.isZeros(candidates[t - p], 0)) {
                CommonOps_DDRM.add(_XList[t - p], candidates[t - p], candidate);
                double l2norm = MyMathUtils.calL2Norm(candidate, TS.getRepairAt(t), dim);
                if (l2norm > this.threshold) {
                    TS.setCandidateAt(t, new DMatrixRMaj(candidate));
                }
            }
        }
    }

    /**
     * Check the monotonicity of a candidate, then use Validate-B() to accept a candidate
     */
    @Override
    protected void evaluate() {
        double l2norm;
        double cosineLowerBound = 0;     // discard a candidate if the cosine is below this value
        int repairIndex;
        Map<Integer, Double> mapRaw = new HashMap<>();
        boolean hasCandidates = false;
        boolean hasRepaired = false;
        // traverse forward candidates
        Iterator<Integer> candIterator = ts.candidateList.iterator();
        while (candIterator.hasNext()) {
            Integer i = candIterator.next();
            if (ts.getCandidateAt(i) == null) {
                throw new RuntimeException("[WARNING]null candidate at " + i);
            }
            double cosine = MyMathUtils.cosineSimilarityOfRepair(ts, i);
            if (cosine < cosineLowerBound) {
                // discard the candidate because it doesn't meet the monotonicity
                candIterator.remove();
                continue;
            }
            if (!ts.getLabelAt(i)) {
                l2norm = MyMathUtils.calL2Norm(ts.getCandidateAt(i), ts.getObsAt(i), dim);
                mapRaw.put(i, l2norm);
            }
        }
        // traverse reverse candidates
        Iterator<Integer> candIteratorR = tsReverse.candidateList.iterator();
        while (candIteratorR.hasNext()) {
            Integer i = candIteratorR.next();
            if (tsReverse.getCandidateAt(i) == null) {
                throw new RuntimeException("[WARNING]null candidate at " + i);
            }
            double cosine = MyMathUtils.cosineSimilarityOfRepair(tsReverse, i);
            if (cosine < cosineLowerBound) {
                candIteratorR.remove();
                continue;
            }
            if (!tsReverse.getLabelAt(i)) {
                l2norm = MyMathUtils.calL2Norm(tsReverse.getCandidateAt(i), tsReverse.getObsAt(i), dim);
                // avoid key conflicts
                mapRaw.put(-(i + 1), l2norm);
            }
        }
        // Sort in descending order by L2 norm and merge the forward and reverse candidates together.
        List<Map.Entry<Integer, Double>> list = new ArrayList<>((mapRaw.entrySet()));
        list.sort(Map.Entry.comparingByValue());
        Map<Integer, Double> map = new LinkedHashMap<>();
        for (Map.Entry<Integer, Double> entry : list) {
            map.put(entry.getKey(), entry.getValue());
        }

        Map<Integer, Double> mapV = new HashMap<>();
        for (Integer integer : map.keySet()) {
            hasCandidates = true;
            repairIndex = integer;
            if (validate(repairIndex, mapV)) {
                if (repairIndex < 0) {
                    repairIndex = -(repairIndex + 1);
                    r = tsLength - repairIndex - 1;
                    tsReverse.setRepairAt(repairIndex, new DMatrixRMaj(tsReverse.getCandidateAt(repairIndex)));
                    ts.setRepairAt(r, new DMatrixRMaj(tsReverse.getCandidateAt(repairIndex)));
                } else {
                    r = repairIndex;
                    int reverseIndex = tsLength - repairIndex - 1;
                    ts.setRepairAt(repairIndex, new DMatrixRMaj(ts.getCandidateAt(repairIndex)));
                    tsReverse.setRepairAt(reverseIndex, new DMatrixRMaj(ts.getCandidateAt(repairIndex)));
                }
                hasRepaired = true;
                break;
            }
        }
        if (!hasRepaired) {
            // There are two scenarios:
            // one is that there are no candidates,
            // and the other is that none of the candidates in the map meet the speed constraints.
            if (hasCandidates) {
                // none of the candidates in the map meet the speed constraints.
                List<Map.Entry<Integer, Double>> list2 = new ArrayList<>((mapV.entrySet()));
                // Sort in descending order of validity.
                list2.sort((o1, o2) -> {
                    return o2.getValue().compareTo(o1.getValue());    // reverse order
                });
                Map<Integer, Double> mapV2 = new LinkedHashMap<>();
                for (Map.Entry<Integer, Double> entry : list2) {
                    mapV2.put(entry.getKey(), entry.getValue());
                }
                Iterator<Integer> iterator2 = mapV2.keySet().iterator();
                Iterator<Double> iteratorValidity = mapV2.values().iterator();
                if (iterator2.hasNext()) {
                    repairIndex = iterator2.next();
                    if (iteratorValidity.next() > 0) {
                        // accept the candidate with the highest positive validity
                        if (repairIndex < 0) {
                            repairIndex = -(repairIndex + 1);
                            r = tsLength - repairIndex - 1;
                            tsReverse.setRepairAt(repairIndex, new DMatrixRMaj(tsReverse.getCandidateAt(repairIndex)));
                            ts.setRepairAt(r, new DMatrixRMaj(tsReverse.getCandidateAt(repairIndex)));
                        } else {
                            r = repairIndex;
                            int reverseIndex = tsLength - repairIndex - 1;
                            ts.setRepairAt(repairIndex, new DMatrixRMaj(ts.getCandidateAt(repairIndex)));
                            tsReverse.setRepairAt(reverseIndex, new DMatrixRMaj(ts.getCandidateAt(repairIndex)));
                        }
                    } else {
                        isConverge = true;
                        System.out.println("Converge with no candidates |= SC");
                    }
                }
                // fixme the else{} part need to be filled
            } else {
                isConverge = true;
                System.out.println("Converge without any candidates.");
            }
        }
        this.RepairedPointSet.add(r);
        ts.clearCandidates();
        tsReverse.clearCandidates();
    }

    /**
     * Check whether a candidate meets the speed constraint s_t and s_g. If it does, accept it; otherwise, provide its validity
     *
     * @param t timesteamp
     * @param m validity Map
     * @return true if meets speed constraint
     */
    @Override
    protected boolean validate(int t, Map<Integer, Double> m) {
        double s_t = Double.MAX_VALUE;
        double s_t_r = Double.MAX_VALUE;
        boolean isReverse = false;
        boolean useBiSC = false;
        double v;
        double vRev = 0;
        int lastLabelP = -1;
        int nextLabelP = -1;
        /*
          t2是对于序列的真实下标，对于正向序列t2=t，对于反向序列，t2=-(t+1);可知t2一定非负
          设置t2的目的在于map中key是不能重复的，反向序列的key要加一取反后存储

          t2 represents the actual index for the sequence.
          For a forward sequence, t2 = t; for a reverse sequence, t2 = -(t + 1).
          It can be inferred that t2 is always non-negative.
          The purpose of setting t2 is to ensure that keys in the map are unique.
          For the reverse sequence, the key is incremented by one and negated before being stored.
         */
        int t2 = t;
        TimeSeries tsUsed = this.ts;
        if (t < 0) {
            // for the reverse candidates
            isReverse = true;
            tsUsed = this.tsReverse;
            t2 = -t - 1;
        }
        if (t2 == 0 || t2 == tsLength - 1) {
            // the first and last data points will not be modified.
            m.put(t, 0.0);
            return false;
        }
        if (isReverse) {
            if (this.RepairedPointSet.contains(tsLength - t2 - 2)) {
                useBiSC = true;
            }
        } else {
            if (this.RepairedPointSet.contains(t2 + 1)) {
                useBiSC = true;
            }
        }
        // Find the two preceding markers and calculate the speed to use as the speed constraint s_t
        for (int i = t2 - 1; i >= 0; i--) {
            if (i == 0) {
                lastLabelP = 0;
                s_t = this.s_g;
                break;
            }
            if (tsUsed.getLabelAt(i)) {
                // Find the nearest preceding markers.
                lastLabelP = i;
                s_t = isReverse ? stRevList[i] : stList[i];
                break;
            }
        }
        if (useBiSC) {
            // Find the two subsequent markers and calculate the speed to use as the speed constraint s_t_r
            for (int i = t2 + 1; i < tsLength; i++) {
                if (i == tsLength - 1) {
                    nextLabelP = tsLength - 1;
                    s_t_r = s_g;
                    break;
                }
                if (tsUsed.getLabelAt(i)) {
                    nextLabelP = i;
                    s_t_r = isReverse ? stList[tsLength - 1 - i] : stRevList[tsLength - 1 - i];
                    break;
                }
            }
            vRev = MyMathUtils.calL2Norm(tsUsed.getCandidateAt(t2), tsUsed.getRepairAt(nextLabelP), dim)
                    / Math.abs((tsUsed.getTimestampAt(t2) - tsUsed.getTimestampAt(nextLabelP)));
        }
        v = MyMathUtils.calL2Norm(tsUsed.getCandidateAt(t2), tsUsed.getRepairAt(lastLabelP), dim)
                / Math.abs((tsUsed.getTimestampAt(t2) - tsUsed.getTimestampAt(lastLabelP)));
        if (useBiSC) {
            // compute validity from both directions
            double minSC = Math.min(Math.min(s_g, s_t), s_t_r);
            if (v > minSC || vRev > minSC) {
                double vo, voRev;
                double validity;
                vo = isReverse ? voRevList[t2] : voList[t2];
                voRev = isReverse ? voList[tsLength - 1 - t2] : voRevList[tsLength - 1 - t2];
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
                long dist = Math.abs(tsUsed.getTimestampAt(t2) - tsUsed.getTimestampAt(lastLabelP));
                long distRev = Math.abs(tsUsed.getTimestampAt(nextLabelP) - tsUsed.getTimestampAt(t2));
                double validityFinal = (validity * Math.exp(distRev) + validityRev * Math.exp(dist)) / (Math.exp(dist) + Math.exp(distRev));
                if (validityFinal > 0) {
                    m.put(t, validityFinal);
                }
                return false;
            }
            return true;
        } else {
            // compute validity from only one direction
            double minSC = Math.min(s_g, s_t);
            if (v > minSC) {
                double vo = isReverse ? voRevList[t2] : voList[t2];
                double validity;
                if (vo <= s_t) {
                    return false;
                }
                double invalidity1 = v - s_t;
                double invalidity2 = vo - s_t;
                validity = Math.max(1 - invalidity1 / invalidity2, 0);
                m.put(t, validity);
                return false;
            }
            return true;
        }
    }

    /**
     * Incrementally update matrix Z and Z^r
     */
    @Override
    protected void updateZ() {
        int row, col, rowBound, bound;
        int colBound = p;
        if (r != -1 && r < tsLength - 1) {
            DMatrixRMaj z = ts.getResidualAt(r);
            int r2 = r + 1;
            if (r2 <= p - 1) {
                // The upper right corner of matrix Z
                rowBound = r2;
                bound = Math.min(rowBound, colBound);
                for (int i = 1; i <= bound; i++) {
                    row = r2 - i + 1;
                    col = (p - i) * dim;
                    for (int j = 0; j < dim; j++) {
                        Z.set(row - 1, col + j, z.get(j, 0));
                    }
                }
            } else {
                // n-p+1<= r2 <= n-1
                rowBound = tsLength - r2;
                bound = Math.min(rowBound, colBound);
                for (int i = 1; i <= bound; i++) {
                    row = r2 - p + i;
                    col = (i - 1) * dim;
                    for (int j = 0; j < dim; j++) {
                        Z.set(row - 1, col + j, z.get(j, 0));
                    }
                }
            }
        }
        // for the reverse series
        int rr = tsLength - 1 - r;
        if (rr != -1 && rr < tsLength - 1) {
            int rr2 = rr + 1;
            DMatrixRMaj zr = tsReverse.getResidualAt(rr);
            if (rr2 <= p - 1) {
                rowBound = rr2;
                bound = Math.min(rowBound, colBound);
                for (int i = 1; i <= bound; i++) {
                    row = rr2 - i + 1;
                    col = (p - i) * dim;
                    for (int j = 0; j < dim; j++) {
                        ZRev.set(row - 1, col + j, zr.get(j, 0));
                    }
                }
            } else {
                rowBound = tsLength - rr2;
                bound = Math.min(rowBound, colBound);
                for (int i = 1; i <= bound; i++) {
                    row = rr2 - p + i;
                    col = (i - 1) * dim;
                    for (int j = 0; j < dim; j++) {
                        ZRev.set(row - 1, col + j, zr.get(j, 0));
                    }
                }
            }
        }
    }
}
