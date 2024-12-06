package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.UniTimeSeries;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.ops.SimpleOperations_DDRM;

import java.util.*;

public class MINOR_BUni extends IMR {
    protected DMatrixRMaj phiReverse;
    protected UniTimeSeries tsReverse;
    protected DMatrixRMaj ZRev;
    protected DMatrixRMaj VRev;
    protected DMatrixRMaj ARev;
    protected DMatrixRMaj BRev;

    protected Set<Integer> RepairedPointSet;

    double f;
    int percentile;
    double sc;

    public MINOR_BUni(UniTimeSeries ts, int p, double threshold, long iterationLimit) {
        this.modelName = "MINOR-B-Uni";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;
        this.iterationLimit = iterationLimit;

        r = -1;
        this.f = 1;
        this.percentile = 99;
        this.tsLength = ts.size();

        Z = new DMatrixRMaj(tsLength - p, p);
        V = new DMatrixRMaj(tsLength - p, 1);
        ZRev = new DMatrixRMaj(tsLength - p, p);
        VRev = new DMatrixRMaj(tsLength - p, 1);
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        ARev = new DMatrixRMaj();
        BRev = new DMatrixRMaj();

        this.phi = new DMatrixRMaj(p, 1);
        this.phiReverse = new DMatrixRMaj(p, 1);

        this.RepairedPointSet = new HashSet<>();
        this.sc = 0;
    }

    @Override
    public RepairResult run() {
        preprocess();
        tsReverse = ts.reverse();
        long startTime = System.nanoTime();
        MatrixPruning();
        int k = 1;
        this.isConverge = false;
        sc = calculatePercentile(calculateSpeeds(), percentile);
        for (int i = 0; i < tsLength; i++) {
            if (this.ts.getLabelAt(i))
                RepairedPointSet.add(i);
        }
        while (true) {
            estimateIC();
            estimateICReverse();
            candidate(this.ts, this.phi);
            candidate(this.tsReverse, this.phiReverse);
            evaluate();
            if (isConverge || k > iterationLimit) {
                System.out.println("[INFO]" + this.modelName + '(' + p + ") stopped at" + k + "th iteration");
                break;
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]time taken: " + duration + "ms");
        double rmse = MyMathUtils.calculateRMSE(ts);
        return new RepairResult(rmse, duration, k);
    }

    /**
     *  repair with Matrix Pruning only
     */
    public RepairResult runMR() {
        preprocess();
        tsReverse = ts.reverse();
        long startTime = System.nanoTime();
        int k = 1;
        this.isConverge = false;
        sc = calculatePercentile(calculateSpeeds(), percentile);
        for (int i = 0; i < tsLength; i++) {
            if (this.ts.getLabelAt(i))
                RepairedPointSet.add(i);
        }
        while (true) {
            MatrixPruning();
            estimateMR();
            estimateMRReverse();
            candidate(this.ts, this.phi);
            candidate(this.tsReverse, this.phiReverse);
            evaluate();
            if (isConverge || k > iterationLimit) {
                System.out.println("[INFO]" + this.modelName + '(' + p + ") stopped at" + k + "th iteration");
                break;
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]time taken: " + duration + "ms");
        double rmse = MyMathUtils.calculateRMSE(ts);
        return new RepairResult(rmse, duration, k);
    }

    @Override
    protected void OLSPre() {
        for (int i = 0; i < tsLength - p; i++) {
            for (int j = 0; j < p; j++) {
                Z.set(i, j, ts.getResidualAt(p - 1 + i - j));
                ZRev.set(i, j, tsReverse.getResidualAt(p - 1 + i - j));
            }
        }
        for (int i = 0; i < tsLength - p; i++) {
            V.set(i, 0, ts.getResidualAt(p + i));
            VRev.set(i, 0, tsReverse.getResidualAt(p + i));
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


    protected void estimateMRReverse() {
        DMatrixRMaj A_inv = new DMatrixRMaj(dim, dim);
        CommonOps_DDRM.invert(ARev, A_inv);
        if (Double.isInfinite(A_inv.get(0, 0)) || Double.isNaN(A_inv.get(0, 0))) {
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(ARev, A_inv);
        }
        CommonOps_DDRM.mult(A_inv, BRev, phiReverse);
    }

    protected void estimateICReverse() {
        if (-1 != r) {
            double subMat, mat1, mat2, mat3, mat4;
            r = tsLength - r - 1;
            r += 1;
            for (int i = 1; i <= p; i++) {
                // a_ii
                if (r >= p + 1 - i && r <= tsLength - i) {
                    subMat = ARev.get(i - 1, i - 1);
                    mat1 = tsReverse.getResidualAt(r - 1) * tsReverse.getResidualAt(r - 1);
                    mat2 = tsReverse.getLastResidualAt(r - 1) * tsReverse.getLastResidualAt(r - 1);
                    mat3 = mat1 - mat2;
                    mat1 = subMat + mat3;
                    ARev.set((i - 1), (i - 1), mat1);
                }
                // a_ij, a_ji, assert i<j
                for (int j = i + 1; j <= p; j++) {
                    if (r >= p + 1 - j && r < p + 1 - i) {
                        subMat = ARev.get(i - 1, j - 1);
                        mat1 = tsReverse.getResidualAt(r - 1) - tsReverse.getLastResidualAt(r - 1);
                        mat2 = tsReverse.getLastResidualAt(r + j - i - 1) * mat1;
                        mat1 = subMat + mat2;
                        ARev.set(j - 1, i - 1, mat1);
                        ARev.set(i - 1, j - 1, mat1);
                    } else if (r > tsLength - j && r <= tsLength - i) {
                        subMat = ARev.get(i - 1, j - 1);
                        mat1 = tsReverse.getResidualAt(r - 1) - tsReverse.getLastResidualAt(r - 1);
                        mat2 = tsReverse.getLastResidualAt(r - j + i - 1) * mat1;
                        mat1 = subMat + mat2;
                        ARev.set(j - 1, i - 1, mat1);
                        ARev.set(i - 1, j - 1, mat1);
                    } else if (r >= p + 1 - i && r <= tsLength - j) {
                        subMat = A.get(i - 1, j - 1);
                        mat1 = tsReverse.getResidualAt(r - 1) - tsReverse.getLastResidualAt(r - 1);
                        mat2 = tsReverse.getLastResidualAt(r + j - i - 1) * mat1;
                        mat3 = tsReverse.getLastResidualAt(r - j + 1 - 1) * mat1;
                        mat4 = mat2 + mat3;
                        mat1 = subMat + mat4;
                        ARev.set(j - 1, i - 1, mat1);
                        ARev.set(i - 1, j - 1, mat1);
                    } // else r<p+1-j || r>n-i, no change
                }
                // bi
                if (r >= p + 1 - i && r < p + 1) {
                    subMat = BRev.get(i - 1, 0);
                    mat1 = tsReverse.getResidualAt(r - 1) - tsReverse.getLastResidualAt(r - 1);
                    mat2 = mat1 * tsReverse.getLastResidualAt(r + i - 1);
                    mat1 = subMat + mat2;
                    BRev.set(i - 1, 0, mat1);
                } else if (r > tsLength - i) {
                    subMat = BRev.get(i - 1, 0);
                    mat1 = tsReverse.getResidualAt(r - 1) - tsReverse.getLastResidualAt(r - 1);
                    mat2 = mat1 * tsReverse.getLastResidualAt(r - i - 1);
                    mat1 = subMat + mat2;
                    BRev.set(i - 1, 0, mat1);
                } else if (r >= p + 1 && r <= tsLength - i) {
                    subMat = BRev.get(i - 1, 0);
                    mat1 = tsReverse.getResidualAt(r - 1) - tsReverse.getLastResidualAt(r - 1);
                    mat2 = mat1 * tsReverse.getLastResidualAt(r + i - 1);
                    mat3 = tsReverse.getLastResidualAt(r - i - 1) * mat1;
                    mat4 = mat2 + mat3;
                    mat1 = subMat + mat4;
                    BRev.set(i - 1, 0, mat1);
                } // else r>p+1-i, no change
            }
            tsReverse.setLastRepairAt(r - 1, tsReverse.getRepairAt(r - 1));
            r -= 1;
            r = tsLength - r - 1;
        }
        DMatrixRMaj A_inv = new DMatrixRMaj(p, p);
        CommonOps_DDRM.invert(ARev, A_inv);
        if (Double.isInfinite(A_inv.get(0, 0)) || Double.isNaN(A_inv.get(0, 0))) {
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(ARev, A_inv);
        }
        CommonOps_DDRM.mult(A_inv, BRev, phiReverse);
    }

    private void candidate(UniTimeSeries TS, DMatrixRMaj Phi) {
        double temp;
        Boolean lb;
        for (int t = p; t < tsLength; t++) {
            lb = TS.getLabelAt(t);
            if (!lb) {
                temp = 0;
                for (int i = 1; i <= p; i++) {
                    temp += Phi.get(i - 1, 0) * TS.getResidualAt(t - i);
                }
                temp += TS.getObsAt(t);
                double l2norm = Math.abs(temp - TS.getRepairAt(t));
                if (l2norm > this.threshold) {
                    TS.setCandidateAt(t, temp);
                }
            }
        }
    }

    @Override
    protected void evaluate() {
        double l2norm;
        int repairIndex;
        Map<Integer, Double> mapRaw = new HashMap<>();
        boolean hasCandidates = false;
        boolean hasRepaired = false;
        Iterator<Integer> candIterator = ts.candidateList.iterator();
        while (candIterator.hasNext()) {
            Integer i = candIterator.next();
            if (ts.getCandidateAt(i) == Double.MAX_VALUE) {
                throw new RuntimeException("[WARNING]null candidate at " + i);
            }
            if (MyMathUtils.isRoundTrip(ts, i)) {
                candIterator.remove();
                continue;
            }
            if (!ts.getLabelAt(i)) {
                l2norm = Math.abs(ts.getCandidateAt(i) - ts.getObsAt(i));
                mapRaw.put(i, l2norm);
            }
        }
        Iterator<Integer> candIteratorR = tsReverse.candidateList.iterator();
        while (candIteratorR.hasNext()) {
            Integer i = candIteratorR.next();
            if (tsReverse.getCandidateAt(i) == Double.MAX_VALUE) {
                throw new RuntimeException("[WARNING]null candidate at " + i);
            }
            if (MyMathUtils.isRoundTrip(tsReverse, i)) {
                candIteratorR.remove();
                continue;
            }
            if (!ts.getLabelAt(i)) {
                l2norm = Math.abs(tsReverse.getCandidateAt(i) - tsReverse.getObsAt(i));
                mapRaw.put(-(i + 1), l2norm);
            }
        }
        List<Map.Entry<Integer, Double>> list = new ArrayList<>((mapRaw.entrySet()));
        list.sort((o1, o2) -> {
            return o2.getValue().compareTo(o1.getValue());    // reverse order, descending order
        });
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
                    tsReverse.setRepairAt(repairIndex, tsReverse.getCandidateAt(repairIndex));
                    ts.setRepairAt(r, tsReverse.getCandidateAt(repairIndex));
                } else {
                    r = repairIndex;
                    int reverseIndex = tsLength - repairIndex - 1;
                    ts.setRepairAt(repairIndex, ts.getCandidateAt(repairIndex));
                    tsReverse.setRepairAt(reverseIndex, ts.getCandidateAt(repairIndex));
                }
                hasRepaired = true;
                break;
            }
        }
        if (!hasRepaired) {
            if (hasCandidates) {
                List<Map.Entry<Integer, Double>> list2 = new ArrayList<>((mapV.entrySet()));
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
                        if (repairIndex < 0) {
                            repairIndex = -(repairIndex + 1);
                            r = tsLength - repairIndex - 1;
                            tsReverse.setRepairAt(repairIndex, tsReverse.getCandidateAt(repairIndex));
                            ts.setRepairAt(r, tsReverse.getCandidateAt(repairIndex));
                        } else {
                            r = repairIndex;
                            int reverseIndex = tsLength - repairIndex - 1;
                            ts.setRepairAt(repairIndex, ts.getCandidateAt(repairIndex));
                            tsReverse.setRepairAt(reverseIndex, ts.getCandidateAt(repairIndex));
                        }
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
        this.RepairedPointSet.add(r);
        ts.clearCandidates();
        tsReverse.clearCandidates();
    }

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
         */
        int t2 = t;
        UniTimeSeries tsUsed = this.ts;
        if (t < 0) {
            isReverse = true;
            tsUsed = this.tsReverse;
            t2 = -t - 1;
        }
        if (t2 == 0 || t2 == tsLength - 1) {
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
        for (int i = t2 - 1; i >= 0; i--) {
            if (i == 0) {
                if (lastLabelP == -1) {
                    lastLabelP = 0;
                    s_t = this.sc;
                } else {
                    s_t = Math.abs((tsUsed.getRepairAt(lastLabelP) - tsUsed.getRepairAt(i)) / (tsUsed.getTimestampAt(lastLabelP) - tsUsed.getTimestampAt(i)));
                }
                break;
            }
            if (!tsUsed.getLabelAt(i)) {
                if (lastLabelP == -1) {
                    lastLabelP = i;
                } else {
                    s_t = Math.abs((tsUsed.getRepairAt(lastLabelP) - tsUsed.getRepairAt(i)) / (tsUsed.getTimestampAt(lastLabelP) - tsUsed.getTimestampAt(i)));
                    break;
                }
            }
        }
        if (useBiSC) {
            for (int i = t2 + 1; i < tsLength; i++) {
                if (i == tsLength - 1) {
                    if (nextLabelP == -1) {
                        nextLabelP = tsLength - 1;
                        s_t_r = sc;
                    } else {
                        s_t = Math.abs((tsUsed.getRepairAt(nextLabelP) - tsUsed.getRepairAt(i)) / (tsUsed.getTimestampAt(nextLabelP) - tsUsed.getTimestampAt(i)));
                    }
                    break;
                }
                if (!tsUsed.getLabelAt(i)) {
                    if (nextLabelP == -1) {
                        nextLabelP = i;
                    } else {
                        s_t = Math.abs((tsUsed.getRepairAt(nextLabelP) - tsUsed.getRepairAt(i)) / (tsUsed.getTimestampAt(nextLabelP) - tsUsed.getTimestampAt(i)));
                        break;
                    }
                }
            }
            vRev = Math.abs(tsUsed.getCandidateAt(t2) - tsUsed.getRepairAt(nextLabelP)
                    / (tsUsed.getTimestampAt(t2) - tsUsed.getTimestampAt(nextLabelP)));
        }
        v = Math.abs((tsUsed.getCandidateAt(t2) - tsUsed.getRepairAt(lastLabelP))
                / (tsUsed.getTimestampAt(t2) - tsUsed.getTimestampAt(lastLabelP)));

        double minSC = Math.min(Math.min(sc, s_t), s_t_r);
        if (v > minSC || vRev > minSC) {
            double vo, voRev;
            double validity;
            vo = Math.abs((tsUsed.getObsAt(t2) - tsUsed.getRepairAt(lastLabelP))
                    / (tsUsed.getTimestampAt(t2) - tsUsed.getTimestampAt(lastLabelP)));
            if (useBiSC) {
                voRev = Math.abs((tsUsed.getObsAt(t2) - tsUsed.getRepairAt(nextLabelP))
                        / (tsUsed.getTimestampAt(t2) - tsUsed.getTimestampAt(nextLabelP)));
                double validityRev;
                double invalidity1 = v - s_t;
                double invalidity2 = vo - s_t;
                if (invalidity1 < 0) {
                    validity = 1;
                } else {
                    if (invalidity2 < 0) {
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
                double validityFinal = (validity * Math.exp(distRev) + validityRev * Math.exp(dist))
                        / (Math.exp(dist) + Math.exp(distRev));
                if (validityFinal > 0) {
                    m.put(t, validityFinal);
                }
            } else {
                if (vo <= s_t) {
                    return false;
                }
                double invalidity1 = v - s_t;
                double invalidity2 = vo - s_t;
                validity = Math.max(1 - invalidity1 / invalidity2, 0);
                m.put(t, validity);
            }
            return false;
        }
        return true;
    }


    protected double[] calculateSpeeds() {
        double[] speeds = new double[tsLength];
        for (int i = 0; i < tsLength - 1; i++) {
            speeds[i] = (Math.abs((ts.getRepairAt(i + 1) - ts.getRepairAt(i)) / (ts.getTimestampAt(i + 1) - ts.getTimestampAt(i))));
        }
        return speeds;
    }

    protected double calculatePercentile(double[] speeds, double percentile) {
        Percentile p = new Percentile();
        p.setData(speeds);
        return p.evaluate(percentile);
    }
}
