package com.MINOR.core;

import com.MINOR.Utils.MyFileUtils;
import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.UniTimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.ops.SimpleOperations_DDRM;

/**
 * 单维时序数据清洗模型
 */
public abstract class UTSCModel extends BaseCleaningModel {
    protected UniTimeSeries ts;

    public UTSCModel() {
        this.dim = 1;
    }

    @Override
    protected void preprocess() {
        for (int i = 0; i < tsLength; i++) {
            boolean lb = ts.getLabelAt(i);
            if (lb) {
                ts.setRepairAt(i, ts.getTruthAt(i));
                ts.setLastRepairAt(i, ts.getTruthAt(i));
            } else {
                ts.setRepairAt(i, ts.getObsAt(i));
                ts.setLastRepairAt(i, ts.getObsAt(i));
            }
        }
    }

    @Override
    protected void OLS() {
        DMatrixRMaj result_temp1 = new DMatrixRMaj(p, p);
        DMatrixRMaj result_temp2 = new DMatrixRMaj(p, tsLength - p);
        DMatrixRMaj result_invert = new DMatrixRMaj(p, p);
        // Z'*Z
        CommonOps_DDRM.multTransA(Z, Z, result_temp1);
        // (Z'*Z)^(-1)*Z'
        CommonOps_DDRM.invert(result_temp1, result_invert);
        if (Double.isInfinite(result_invert.get(0, 0)) || Double.isNaN(result_invert.get(0, 0))) {
            // for singular matrix (Z'Z), use pseudo-invert instead.
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(result_temp1, result_invert);
        }
        CommonOps_DDRM.multTransB(result_invert, Z, result_temp2);
        CommonOps_DDRM.mult(result_temp2, V, phi);
    }

    @Override
    public double getRMSEDirty() {
        if (this.ts == null) {
            System.out.println("[WARNING]null time series when calculating RMSE of dirty data.");
            return Double.NaN;
        } else {
            return MyMathUtils.calculateRMSE_dirty(this.ts);
        }
    }

    @Override
    public double getRMSE() {
        if (this.ts == null) {
            System.out.println("[WARNING]null time series when calculating RMSE of repaired data.");
            return Double.NaN;
        } else {
            return MyMathUtils.calculateRMSE(this.ts);
        }
    }

    @Override
    public void exportTimeSeries(String fileName) {
        MyFileUtils.UTSExporter(fileName, this.ts);
    }
}
