package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

import java.util.ArrayList;
import java.util.List;

public class VARX extends MTSCModel {
    public VARX(){

    }

    public VARX(TimeSeries ts, int p, double threshold){
        this.modelName = "VARX";

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

    /**
     * used for application experiment
     *
     * @return cleaned data
     */
    public TimeSeries runApp() {
        preprocess();
        OLSPre();
        OLS();
        repair();
        return this.ts;
    }

    @Override
    protected void OLSPre() {
        List<DMatrixRMaj> z = new ArrayList<>();
        for (int i = 0; i < tsLength; i++) {
            DMatrixRMaj resultMatrix = new DMatrixRMaj(dim, 1);
            CommonOps_DDRM.subtract(ts.getRepairAt(i), ts.getObsAt(i), resultMatrix);
            z.add(resultMatrix);
        }
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
    }

    protected void repair() {
        DMatrixRMaj temp = new DMatrixRMaj(dim, 1);
        DMatrixRMaj phiSubMat = new DMatrixRMaj(dim, dim);
        DMatrixRMaj z = new DMatrixRMaj(dim, 1);
        DMatrixRMaj z2 = new DMatrixRMaj(dim, 1);
        for (int t = p; t < tsLength; t++) {
            if (!ts.getLabelAt(t)) {
                temp.zero();
                for (int i = 1; i <= p; i++) {
                    // Obtain the corresponding submatrix when the lag coefficient is p
                    CommonOps_DDRM.extract(phi, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, phiSubMat);
                    CommonOps_DDRM.subtract(ts.getRepairAt(t - i), ts.getObsAt(t - i), z);
                    CommonOps_DDRM.multTransA(phiSubMat, z, z2);
                    CommonOps_DDRM.add(temp, z2, temp);
                }
                CommonOps_DDRM.add(temp, ts.getObsAt(t), temp);
                // The L2 norm between the predicted value and the original predicted value is greater than the threshold
                double l2norm = MyMathUtils.calL2Norm(temp, ts.getRepairAt(t), dim);
                if (l2norm > this.threshold) {
                    ts.setRepairAt(t, new DMatrixRMaj(temp));
                }
            }
        }
    }
}
