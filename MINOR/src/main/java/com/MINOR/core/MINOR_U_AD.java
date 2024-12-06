package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;

@Deprecated
public class MINOR_U_AD extends MINOR_U implements AnomalyDetector{
    public MINOR_U_AD(TimeSeries ts, int p, double threshold, long iterationLimit) {
        super(ts,p,threshold,iterationLimit);
        this.modelName = "MINOR_U_AD";
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        preprocess();
        MatrixPruning();
        int k = 1;
        this.isConverge = false;
        s_g = calculatePercentile(calculateSpeeds(this.f), this.percentile);
        this.detectAnomalies(ts);
        this.ADEvaluate(ts);
        while (true) {
            estimateIC();
            candidate();
            evaluate();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + this.modelName + "_pvad(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        ts.revertAD();
        return new RepairResult(MyMathUtils.calculateRMSE(ts), duration, k);
    }
}
