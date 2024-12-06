package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.SortedList;
import com.MINOR.entity.TimePoint;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.ops.SimpleOperations_DDRM;

import java.util.Comparator;

public class MINOR_O extends MINOR_U {
    SortedList<Double> speedList;
    double s_t;     // The speed of the last two points in the corrected sequence.

    TimeSeries tsTotal;
    int tsLengthTotal;
    int start;

    public MINOR_O(TimeSeries _ts) {
        this.tsTotal = _ts;
        this.threshold = 0;
        this.p = 1;
        this.tsLengthTotal = tsTotal.size();
        this.dim = _ts.getDimension();
        r = start - 1;
        s_t = 0;
        this.modelName = "MINOR-O";

        Z = new DMatrixRMaj(tsLengthTotal - p, dim * p);
        V = new DMatrixRMaj(tsLengthTotal - p, dim);
        this.phi = new DMatrixRMaj(p * dim, dim);

        speedList = new SortedList<>(Comparator.reverseOrder());

        start = 3;
        this.ts = new TimeSeries();
        this.tsLength = 0;
    }

    @Override
    public RepairResult run() {
        preprocess();
        long startTime = System.nanoTime();
        OLSPre();
        repair();
        long endTime = System.nanoTime();
        System.out.println("[INFO]" + modelName + "(" + 1 + ") Done.");
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        double rmse = MyMathUtils.calculateRMSE(ts);
        return new RepairResult(rmse, duration, 1);
    }

    @Override
    protected void preprocess() {
        // assign the first 3 points to it
        for (int i = 0; i < start; i++) {
            TimePoint tp = tsTotal.getP(i);
            tp.setValRepaired(tp.getValTruth());
            tp.setValLastRepaired(tp.getValTruth());
            this.ts.add(tp);
        }
        tsLength += start;
        A = new DMatrixRMaj(p * dim, p * dim);
        B = new DMatrixRMaj(p * dim, dim);
    }

    protected void repair() {
        for (int i = start; i < tsLengthTotal; i++) {
            TimePoint tp = tsTotal.getP(i);
            if (!tp.getLabel()) {
                s_t = MyMathUtils.calSpeed(ts.getP(tsLength - 1), ts.getP(tsLength - 2));
                double v, vo, deltaT;
                DMatrixRMaj z = new DMatrixRMaj(dim, 1);
                DMatrixRMaj modify = new DMatrixRMaj(dim, 1);
                DMatrixRMaj yhat = new DMatrixRMaj(dim, 1);
                CommonOps_DDRM.subtract(ts.getRepairAt(i - 1), ts.getObsAt(i - 1), z);
                CommonOps_DDRM.multTransA(phi, z, modify);
                CommonOps_DDRM.add(modify, tp.getValObs(), yhat);
                deltaT = Math.abs(tp.getTimestamp() - ts.getTimestampAt(tsLength - 1));
                v = MyMathUtils.calL2Norm(yhat, ts.getRepairAt(tsLength - 1), this.dim)
                        / deltaT;
                if (v <= s_t) {
                    tp.setValRepaired(yhat);
                } else {
                    vo = MyMathUtils.calL2Norm(tp.getValObs(), ts.getRepairAt(tsLength - 1), this.dim)
                            / deltaT;
                    if (v < vo) {
                        tp.setValRepaired(yhat);
                    } else {
                        tp.setValRepaired(tp.getValObs());
                    }
                }
                v = MyMathUtils.calL2Norm(tp.getValRepaired(), ts.getRepairAt(tsLength - 1), this.dim)
                        / deltaT;
            } else {
                tp.setValRepaired(tp.getValTruth());
                double v = MyMathUtils.calL2Norm(tp.getValTruth(), ts.getRepairAt(tsLength - 1), this.dim)
                        / Math.abs(tp.getTimestamp() - ts.getTimestampAt(tsLength - 1));
            }
            tp.setValLastRepaired(tp.getValObs());
            ts.add(new TimePoint(tp));
            tsLength++;
            r = i;
            estimateIC();
        }
    }

    @Override
    protected void estimateIC() {
        if (-1 != r) {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    double a = A.get(i, j);
                    double b = B.get(i, j);
                    A.set(i, j, a + ts.getResidualAt(r).get(i, 0) * ts.getResidualAt(r).get(j, 0));
                    B.set(i, j, b + ts.getResidualAt(r - 1).get(i, 0) * ts.getResidualAt(r - 1).get(j, 0));
                }
            }
        }
        DMatrixRMaj A_inv = new DMatrixRMaj(dim * p, dim * p);
        CommonOps_DDRM.invert(A, A_inv);
        if (Double.isInfinite(A_inv.get(0, 0)) || Double.isNaN(A_inv.get(0, 0))) {
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(A, A_inv);
        }
        CommonOps_DDRM.mult(A_inv, B, phi);
    }
}
