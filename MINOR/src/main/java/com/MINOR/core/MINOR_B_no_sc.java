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

public class MINOR_B_no_sc extends MINOR_B {
    public MINOR_B_no_sc() {

    }

    public MINOR_B_no_sc(TimeSeries ts, int p, double threshold, long iterationLimit) {
        this.modelName = "MINOR-B-w/oSC";

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
     * Check the monotonicity of a candidate, then accept a candidate
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

        for (Integer integer : map.keySet()) {
            repairIndex = integer;
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
            break;
        }
        this.RepairedPointSet.add(r);
        ts.clearCandidates();
        tsReverse.clearCandidates();
    }
}