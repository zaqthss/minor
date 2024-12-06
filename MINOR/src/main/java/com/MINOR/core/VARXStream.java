package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimePoint;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

public class VARXStream extends VARX {
    TimeSeries tsTotal;
    int tsLengthTotal;
    int start;

    public VARXStream(TimeSeries _ts, double threshold) {
        this.modelName = "VARXStream";

        this.tsTotal = _ts;
        this.threshold = threshold;
        this.p = 1;
        this.tsLengthTotal = tsTotal.size();
        this.dim = _ts.getDimension();

        Z = new DMatrixRMaj(tsLengthTotal - p, dim * p);
        V = new DMatrixRMaj(tsLengthTotal - p, dim);
        this.phi = new DMatrixRMaj(p * dim, dim);

        start = 3;
    }

    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocess();
        OLSPre();
        repair();
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        System.out.println("[INFO]" + this.modelName + '(' + 1 + ") Done.");
        return new RepairResult(MyMathUtils.calculateRMSE(ts), duration, 1);
    }

    @Override
    protected void preprocess() {
        this.ts = new TimeSeries();
        // assign the first 3 points to it
        for (int i = 0; i < start; i++) {
            TimePoint tp = tsTotal.getP(i);
            tp.setValRepaired(tp.getValTruth());
            tp.setValLastRepaired(tp.getValTruth());
            this.ts.add(tp);
        }
        this.tsLength = 3;
    }

    protected void updateZV(int r) {
        if (r > tsLengthTotal - 2) {
            return;
        }
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < dim; k++) {
                Z.set(r, dim * j + k, ts.getResidualAt(r).get(k, 0));
                V.set(r - 1, k, ts.getResidualAt(r - 1).get(k, 0));
            }
        }
    }

    @Override
    protected void repair() {
        for (int i = start; i < tsLengthTotal; i++) {
            TimePoint tp = tsTotal.getP(i);
            if (!tp.getLabel()) {
                OLS();
                DMatrixRMaj z = new DMatrixRMaj(dim, 1);
                DMatrixRMaj modify = new DMatrixRMaj(dim, 1);
                DMatrixRMaj yhat = new DMatrixRMaj(dim, 1);
                CommonOps_DDRM.subtract(ts.getRepairAt(i - 1), ts.getObsAt(i - 1), z);
                CommonOps_DDRM.multTransA(phi, z, modify);
                CommonOps_DDRM.add(modify, tp.getValObs(), yhat);
                double l2norm = MyMathUtils.calL2Norm(yhat, tp.getValObs(), dim);
                if (l2norm > this.threshold) {
                    tp.setValRepaired(yhat);
                } else {
                    tp.setValRepaired(tp.getValObs());
                }
            } else {
                tp.setValRepaired(tp.getValTruth());
            }
            ts.add(new TimePoint(tp));
            tsLength++;
            updateZV(i);
        }
    }
}
