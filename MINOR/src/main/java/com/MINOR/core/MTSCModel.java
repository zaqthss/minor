package com.MINOR.core;

import com.MINOR.Utils.MyFileUtils;
import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.ops.SimpleOperations_DDRM;

/**
 * 多维时序数据清洗模型
 */
public abstract class MTSCModel extends BaseCleaningModel {
    protected TimeSeries ts;

    MTSCModel() {
    }

    /**
     * Replace the marked data points' observed values with true values, and restore values such as lastRepair and lastModify,
     * so that different settings can be tested consecutively during testing.
     */
    @Override
    protected void preprocess() {
        for (int i = 0; i < tsLength; i++) {
            DMatrixRMaj mat;
            DMatrixRMaj obs = ts.getObsAt(i);
            DMatrixRMaj truth = ts.getTruthAt(i);
            if (ts.getLabelAt(i)) {
                mat = new DMatrixRMaj(truth);
            } else {
                mat = new DMatrixRMaj(obs);
            }
            ts.setRepairAt(i, new DMatrixRMaj(mat));
            ts.setLastRepairAt(i, new DMatrixRMaj(mat));
            ts.setLastModifyAt(i, new DMatrixRMaj(mat));
        }
    }

    /**
     * Calculate the least squares solution using the Moore-Penrose pseudoinverse matrix to obtain the model parameters phi,
     * with all vectors in column vector form.
     */
    @Override
    protected void OLS() {
        DMatrixRMaj result_temp1 = new DMatrixRMaj(dim * p, dim * p);
        DMatrixRMaj result_temp2 = new DMatrixRMaj(dim * p, tsLength - p);
        DMatrixRMaj result_invert = new DMatrixRMaj(dim * p, dim * p);
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
        MyFileUtils.MTSExporter(fileName, this.ts);
    }
}
