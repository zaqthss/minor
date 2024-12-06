package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

public class VAR extends MTSCModel {
    public VAR(TimeSeries ts, int p, double threshold) {
        this.modelName = "VAR";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;

        this.dim = ts.getDimension();
        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, dim * p);
        V = new DMatrixRMaj(tsLength - p, dim);
        this.phi = new DMatrixRMaj(p * dim, dim);
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
                for (int k = 0; k < dim; k++) {
                    Z.set(i, dim * j + k, ts.getRepairAt(p - 1 + i - j).get(k, 0));
                }
            }
        }
        for (int i = 0; i < tsLength - p; i++) {
            for (int k = 0; k < dim; k++) {
                V.set(i, k, ts.getRepairAt(p + i).get(k, 0));
            }
        }
    }

    private void repair() {
        DMatrixRMaj temp = new DMatrixRMaj(dim, 1);
        DMatrixRMaj phiSubMat = new DMatrixRMaj(dim, dim);
        DMatrixRMaj z = new DMatrixRMaj(dim, 1);
        for (int t = p; t < tsLength; t++) {
            if (!ts.getLabelAt(t)) {
                temp.zero();
                for (int i = 1; i <= p; i++) {
                    // Obtain the corresponding submatrix when the lag coefficient is p
                    CommonOps_DDRM.extract(phi, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, phiSubMat);
                    CommonOps_DDRM.multTransA(phiSubMat, ts.getRepairAt(t - i), z);
                    CommonOps_DDRM.add(temp, z, temp);
                }
                // The L2 norm between the predicted value and the original predicted value is greater than the threshold
                double l2norm = MyMathUtils.calL2Norm(temp, ts.getRepairAt(t), dim);
                if (l2norm > this.threshold) {
                    ts.setRepairAt(t, new DMatrixRMaj(temp));
                }
            }
        }
    }
}
