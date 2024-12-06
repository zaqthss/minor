package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;

@Deprecated
public class MINOR_B_AD extends MINOR_B implements AnomalyDetector {
    public MINOR_B_AD(TimeSeries ts, int p, double threshold, long iterationLimit) {
        super(ts, p, threshold, iterationLimit);
        this.modelName = "MINOR_B_AD";
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocess();
        tsReverse = ts.reverse();
        MatrixPruning();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        for (int i = 0; i < tsLength; i++) {
            if (ts.getLabelAt(i)) {
                RepairedPointSet.add(i);
            }
        }
        this.detectAnomalies(ts);
        this.ADEvaluate(ts);
        while (true) {
            estimateIC();
            estimateICReverse();
            candidate(this.ts, this.phi);
            candidate(this.tsReverse, this.phiReverse);
            evaluate();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + this.modelName + "_bipvad(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            if (k % (iterationLimit / 10) == 0) {
                System.out.println(k);
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        ts.revertAD();
        return new RepairResult(MyMathUtils.calculateRMSE(ts), duration, k);
    }
}
