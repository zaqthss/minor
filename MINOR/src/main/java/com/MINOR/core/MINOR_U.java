package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.MatrixFeatures_DDRM;

import java.util.*;

public class MINOR_U extends MINOR_base {
    protected double s_g;   // global speed constraint
    protected double percentile;
    protected double f;
    protected int fistLabel;
    protected double[] stList; // 每个标记点与前一个标记点的速度
    protected double[] voList;    // 每个点观测值与前一个点修复值速度
    protected DMatrixRMaj X;
    protected DMatrixRMaj[] XList;

    public MINOR_U(TimeSeries ts, int p, double threshold, long iterationLimit) {
        this.modelName = "MINOR-U";

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
        this.phi = new DMatrixRMaj(p * dim, dim);

        this.s_g = 0;
        this.stList = new double[tsLength];
        this.voList = new double[tsLength];
    }

    public MINOR_U() {
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocessIC();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        while (true) {
            estimateIC();
            candidateIC();
            evaluate();
            updateZ();
            if (isConverge || k >= iterationLimit) {
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

    /**
     * used for application experiment
     *
     * @return cleaned data
     */
    public TimeSeries runApp() {
        preprocessIC();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        while (true) {
            estimateIC();
            candidateIC();
            evaluate();
            updateZ();
            if (isConverge || k >= iterationLimit) {
                break;
            }
            k++;
        }
        return ts;
    }

    /**
     *  print time cost of each step when running
     */
    public RepairResult runWithTimer() {
        double[] timeCost = new double[4];
        long startTime = System.nanoTime();
        preprocessIC();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        while (true) {
            long startTime1 = System.nanoTime();
            estimateIC();
            long endTime1 = System.nanoTime();
            candidateIC();
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
        System.out.println("[INFO]time cost of 4 parts: estimate: " + timeCost[0]
                + ", candidate: " + timeCost[1]
                + ", evaluate: " + timeCost[2]
                + ", updateZ: " + timeCost[3]);
        return result;
    }

    /**
     *  repair without any Incremental Computation
     */
    public RepairResult runSD() {
        long startTime = System.nanoTime();
        preprocessIC();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        while (true) {
            estimate();
            candidate();
            evaluate();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + this.modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            if (k % (iterationLimit / 10) == 0) {
                System.out.println("[REPAIRING] " + k + " / " + iterationLimit);
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
     *  repair without Incremental Candidate Generation
     */
    public RepairResult runWithoutCandiIC() {
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        long startTime = System.nanoTime();
        preprocessIC();
        MatrixPruning();
//        preprocessIC();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        double[] timeCost = new double[3];
        while (true) {
            long startTime1 = System.nanoTime();
            estimateIC();
            long endTime1 = System.nanoTime();
            candidate();
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
        System.out.println("[INFO]time cost of 4 parts: estimate: " + timeCost[0]
                + ", candidate: " + timeCost[1]
                + ", evaluate: " + timeCost[2]);
        return result;
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
//                voList[i] = CustomMathUtils.calL2Norm(obs, ts.getRepairAt(i - 1), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(i - 1));
                if (lastLabel == -1) {
                    voList[i] = MyMathUtils.calL2Norm(obs, ts.getRepairAt(0), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(0));
                } else {
                    voList[i] = MyMathUtils.calL2Norm(obs, ts.getTruthAt(lastLabel), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(lastLabel));
                }
            }
        }
    }

    /**
     * initialize data, memorize vo and st,
     * and initialize variables for Incremental Candidate Generation
     */
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
                    // The speed between labeled data points will not change after preprocessing.
                    stList[i] = MyMathUtils.calL2Norm(truth, ts.getTruthAt(lastLabel), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(lastLabel));
                } else {
                    this.fistLabel = i;
                    stList[i] = MyMathUtils.calL2Norm(truth, ts.getRepairAt(0), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(0));
                }
                lastLabel = i;
            }
            if (i > 0) {
//                voList[i] = CustomMathUtils.calL2Norm(obs, ts.getRepairAt(i - 1), dim) / Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(i - 1));
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
     * generate candidates with Incremental Computation of Z
     */
    protected void candidateIC() {
        // Y_hat = X + Phi^T * Z^T
        // each column of Y_hat is a candidate at the index
        DMatrixRMaj[] candidates = new DMatrixRMaj[tsLength - p];
        DMatrixRMaj candidate = new DMatrixRMaj(dim, 1);
        CommonOps_DDRM.columnsToVector(CommonOps_DDRM.multTransAB(phi, Z, null), candidates);
        for (int t = p; t < tsLength; t++) {
            if (!ts.getLabelAt(t) && !MatrixFeatures_DDRM.isZeros(candidates[t - p], 0)) {
                CommonOps_DDRM.add(XList[t - p], candidates[t - p], candidate);
                double l2norm = MyMathUtils.calL2Norm(candidate, ts.getRepairAt(t), dim);
                if (l2norm > this.threshold) {
                    ts.setCandidateAt(t, new DMatrixRMaj(candidate));
                }
            }
        }
    }

//    @Override
//    protected void candidate(){
//        OLSPre();
//        DMatrixRMaj[] candidates = new DMatrixRMaj[tsLength - p];
//        DMatrixRMaj candidate = new DMatrixRMaj(dim, 1);
//        CommonOps_DDRM.columnsToVector(CommonOps_DDRM.multTransAB(phi, Z, null), candidates);
//        for (int t = p; t < tsLength; t++) {
//            if (!ts.getLabelAt(t)) {
//                CommonOps_DDRM.add(XList[t - p], candidates[t - p], candidate);
//                double l2norm = MyMathUtils.calL2Norm(candidate, ts.getRepairAt(t), dim);
//                if (l2norm > this.threshold) {
//                    ts.setCandidateAt(t, new DMatrixRMaj(candidate));
//                }
//            }
//        }
//    }

    /**
     * use Validate-U() to accept a candidate
     */
    @Override
    protected void evaluate() {
        double l2norm;
        int repairIndex;
        Map<Integer, Double> mapRaw = new HashMap<>();
        boolean hasCandidates = false;
        boolean hasRepaired = false;
        for (Integer i : ts.candidateList) {
            if (ts.getCandidateAt(i) == null) {
                throw new RuntimeException("[WARNING]null candidate at " + i);
            }
            if (!ts.getLabelAt(i)) {
                l2norm = MyMathUtils.calL2Norm(ts.getCandidateAt(i), ts.getObsAt(i), dim);
                mapRaw.put(i, l2norm);
            }
        }
        // Sort in ascending order based on L2 norm
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
                ts.setRepairAt(repairIndex, new DMatrixRMaj(ts.getCandidateAt(repairIndex)));
                r = repairIndex;
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
                list2.sort((o1, o2) -> o2.getValue().compareTo(o1.getValue()));
                Map<Integer, Double> mapV2 = new LinkedHashMap<>();
                for (Map.Entry<Integer, Double> entry : list2) {
                    mapV2.put(entry.getKey(), entry.getValue());
                }
                Iterator<Integer> iterator2 = mapV2.keySet().iterator();
                Iterator<Double> iteratorValidity = mapV2.values().iterator();
                if (iterator2.hasNext()) {
                    // accept the candidate with the highest positive validity
                    repairIndex = iterator2.next();
                    if (iteratorValidity.next() > 0) {
                        ts.setRepairAt(repairIndex, new DMatrixRMaj(ts.getCandidateAt(repairIndex)));
                        if (r != -1) {
                            ts.setLastRepairAt(r, new DMatrixRMaj(ts.getRepairAt(r)));
                        }
                        r = repairIndex;
                    } else {
                        isConverge = true;
                        System.out.println("Converge with no candidates |= SC");
                    }
                }
            } else {
                isConverge = true;
                System.out.println("Converge without any candidates.");
            }
        }
        ts.clearCandidates();
    }

    /**
     * Check whether a candidate meets the speed constraint s_t and s_g. If it does, accept it; otherwise, provide its validity
     *
     * @param t timesteamp
     * @param m validity Map
     * @return true if meets speed constraint
     */
    protected boolean validate(int t, Map<Integer, Double> m) {
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
                return false;
            } else {
                double invalidity1 = v - s_t;
                double invalidity2 = vo - s_t;
                double validity = Math.max(1 - invalidity1 / invalidity2, 0);
                m.put(t, validity);
            }
            return false;
        }
        return true;
    }

    /**
     * Incrementally update matrix Z
     */
    protected void updateZ() {
        if (r != -1 && r < tsLength - 1) {
            if (p > 1 && r < p - 1) {
                return;
            }
            DMatrixRMaj z = ts.getResidualAt(r);
            int r2 = r + 1;
            int row, col;
            int rowBound = tsLength - r2;
            int colBound = p;
            int bound = Math.min(rowBound, colBound);
            for (int i = 1; i <= bound; i++) {
                row = r2 - p + i;
                col = (i - 1) * dim;
                for (int j = 0; j < dim; j++) {
                    Z.set(row - 1, col + j, z.get(j, 0));
                }
            }
        }
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
