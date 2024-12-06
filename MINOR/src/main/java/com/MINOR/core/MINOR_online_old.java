package com.MINOR.core;

import com.MINOR.Utils.MyFileUtils;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

@Deprecated
public class MINOR_online_old extends MINOR_base {
    public MINOR_online_old(int dim) {
        this.dim = dim;
        this.ts = null;
        this.p = 1;
        this.phi = new DMatrixRMaj(dim, dim);
        r = -1;
    }

    public TimeSeries run(TimeSeries subTs) {
        preprocess(subTs);
        if (this.ts == null) {
            // 需要保证第一段子序列全为标记数据
            this.ts = subTs;
            tsLength = subTs.size();
            if (tsLength >= 2) {
                Z = new DMatrixRMaj(tsLength - p, dim * p);
                V = new DMatrixRMaj(tsLength - p, dim);
                A = new DMatrixRMaj();
                B = new DMatrixRMaj();
                MatrixPruning();
                estimateMR();
            }
        } else if (tsLength < 2) {
            this.ts.concat(subTs);
            tsLength += subTs.size();
            Z = new DMatrixRMaj(tsLength - p, dim * p);
            V = new DMatrixRMaj(tsLength - p, dim);
            A = new DMatrixRMaj();
            B = new DMatrixRMaj();
            MatrixPruning();
            estimateMR();
        } else {
            repair(subTs);
            this.ts.concat(subTs);
            tsLength += subTs.size();
            Z.reshape(tsLength - p, dim * p);
            V.reshape(tsLength - p, dim);
            MatrixPruning();
            estimateMR();
        }
        // 获得phi，本轮的phi将提供给下一轮repair
        return this.ts;
    }

    public void printResult(String resultFilename, String phiFilename) {
        MyFileUtils.MTSExporter(resultFilename, ts);
        MyFileUtils.PhiExporter(phiFilename, this.phi);
    }

    /**
     * 给时序数据片段中数据点的ValRepaired属性赋值
     *
     * @param subTs 到来的时序数据片段
     */
    private void preprocess(TimeSeries subTs) {
        int n = subTs.size();
        for (int i = 0; i < n; i++) {
            DMatrixRMaj mat = new DMatrixRMaj(dim, 1);
            DMatrixRMaj obs = subTs.getObsAt(i);
            DMatrixRMaj truth = subTs.getTruthAt(i);
            if (ts == null || ts.getLabelAt(i)) {
//                mat.set(dim, 0, truth.get(0, 0));
                mat = new DMatrixRMaj(truth);
            } else {
//                mat.set(dim, 0, obs.get(0, 0));
                mat = new DMatrixRMaj(obs);
            }
            subTs.setRepairAt(i, new DMatrixRMaj(mat));
            subTs.setLastRepairAt(i, new DMatrixRMaj(mat));
        }
    }


    private void repair(TimeSeries subTs) {
        int n = subTs.size();
        DMatrixRMaj temp = new DMatrixRMaj(dim, 1);
        DMatrixRMaj phiPowMat = new DMatrixRMaj(phi);
        DMatrixRMaj phiPowMatTemp = new DMatrixRMaj(dim, dim);
        DMatrixRMaj z;
        DMatrixRMaj z2 = new DMatrixRMaj(dim, 1);
        for (int t = 0; t < n; t++) {
            if (!ts.getLabelAt(t)) {
                temp.zero();
                for (int i = n - t - 1; i < n - 1; i++) {
                    CommonOps_DDRM.mult(phiPowMat, phi, phiPowMatTemp);
                    phiPowMat = new DMatrixRMaj(phiPowMatTemp);
                }
                z = ts.getResidualAt(tsLength - 1);
                CommonOps_DDRM.multTransA(phiPowMat, z, z2);
                CommonOps_DDRM.add(temp, z2, temp);
                CommonOps_DDRM.add(temp, subTs.getObsAt(t), temp);
                subTs.setRepairAt(t, new DMatrixRMaj(temp));
                r = t;
            }
        }
    }
}
