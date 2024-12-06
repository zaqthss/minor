package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.UniTimeSeries;
import org.ejml.data.DMatrixRMaj;

public class AR extends UTSCModel {
    public AR(UniTimeSeries ts, int p, double threshold) {
        this.modelName = "AR";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;

        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, p);
        V = new DMatrixRMaj(tsLength - p, 1);
        this.phi = new DMatrixRMaj(p, 1);
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocess();
        OLSPre();
        OLS();
        repair();
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]" + this.modelName + "(" + p + ") Done.");
        return new RepairResult(MyMathUtils.calculateRMSE(ts), duration, 1);
    }

    @Override
    protected void OLSPre() {
        for (int i = 0; i < tsLength - p; i++) {
            for (int j = 0; j < p; j++) {
                Z.set(i, j, ts.getRepairAt(p - 1 + i - j));
            }
        }
        for (int i = 0; i < tsLength - p; i++) {
            V.set(i, 0, ts.getRepairAt(p + i));
        }
    }

    private void repair() {
        double temp;
        Boolean lb;
        for (int t = p; t < tsLength; t++) {
            lb = ts.getLabelAt(t);
            if (!lb) {
                temp = 0;
                for (int i = 1; i <= p; i++) {
                    temp += phi.get(i - 1, 0) * ts.getRepairAt(t - i);
                }
                double l2norm = Math.abs(temp - ts.getRepairAt(t));
                if (l2norm > this.threshold) {
                    ts.setRepairAt(t, temp);
                }
            }
        }
    }
}
