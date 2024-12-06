package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.ops.SimpleOperations_DDRM;

import java.util.ArrayList;
import java.util.List;

/**
 * 基础的修复方法, 没有anomaly detection也没有post validation
 */
public class MINOR_base extends MTSCModel {
    protected long iterationLimit;

    protected DMatrixRMaj A;
    protected DMatrixRMaj B;
    protected boolean isConverge;
    protected int r;

    public MINOR_base(TimeSeries ts, int p, double threshold, long iterationLimit) {
        this.modelName = "MINOR-base";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;
        this.iterationLimit = iterationLimit;

        r = -1;

        this.isConverge = false;
        this.dim = ts.getDimension();
        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, dim * p);
        V = new DMatrixRMaj(tsLength - p, dim);
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        this.phi = new DMatrixRMaj(p * dim, dim);
    }

    public MINOR_base() {
    }


    @Override
    public RepairResult run() {
        long startTime = System.nanoTime();
        preprocess();
        MatrixPruning();
        int k = 1;
        this.isConverge = false;
        while (true) {
            estimateIC();
            candidate();
            evaluate();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done. " + modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        double rmse = MyMathUtils.calculateRMSE(ts);
        return new RepairResult(rmse, duration, k);
    }

    /**
     *  repair with Matrix Pruning only
     */
    public RepairResult runMR() {
        long startTime = System.nanoTime();
        preprocess();
        int k = 1;
        this.isConverge = false;
        while (true) {
            MatrixPruning();
            estimateMR();
            candidate();
            evaluate();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        double rmse = MyMathUtils.calculateRMSE(ts);
        return new RepairResult(rmse, duration, k);
    }

    /**
     * repair without Matrix Pruning and Incremental Computation
     */
    public RepairResult runSD() {
        long startTime = System.nanoTime();
        preprocess();
        int k = 1;
        this.isConverge = false;
        while (true) {
            estimate();
            candidate();
            evaluate();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + modelName + "(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            k++;
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        double rmse = MyMathUtils.calculateRMSE(ts);
        return new RepairResult(rmse, duration, k);
    }

    @Override
    protected void OLSPre() {
        List<DMatrixRMaj> z = new ArrayList<>();
        for (int i = 0; i < tsLength; i++) {
            DMatrixRMaj resultMatrix = new DMatrixRMaj(dim, 1);
            CommonOps_DDRM.subtract(ts.getRepairAt(i), ts.getObsAt(i), resultMatrix);
            z.add(resultMatrix);
        }
        // Matrix Z
        for (int i = 0; i < tsLength - p; i++) {
            for (int j = 0; j < p; j++) {
                for (int k = 0; k < dim; k++) {
                    Z.set(i, dim * j + k, z.get(p - 1 + i - j).get(k, 0));
                }
            }
        }
        // Matrix V
        for (int i = 0; i < tsLength - p; i++) {
            for (int k = 0; k < dim; k++) {
                V.set(i, k, z.get(p + i).get(k, 0));
            }
        }
    }

    /**
     * remove rows with all zeros of Z and the corresponding rows of V, and generate matrix A and B for OLS
     */
    protected void MatrixPruning() {
        OLSPre();
        DMatrixRMaj[] wrapper = {Z, V};
        MyMathUtils.MatrixPruning(wrapper);
//        Z = wrapper[0];
//        V = wrapper[1];
        A.reshape(Z.numCols, Z.numCols);
        B.reshape(Z.numCols, V.numCols);
        CommonOps_DDRM.multTransA(wrapper[0], wrapper[0], A);         // A_0 = Z'*Z
        CommonOps_DDRM.multTransA(wrapper[0], wrapper[1], B);         // B_0 = Z'*V
    }

    /**
     * learn parameter Phi without Matrix Pruning and Incremental Computation
     */
    protected void estimate() {
        OLSPre();
        OLS();
    }

    /**
     * learn parameter Phi with Matrix Pruning only
     */
    protected void estimateMR() {
        DMatrixRMaj A_inv = new DMatrixRMaj(dim * p, dim * p);
        // (A)^(-1)
        CommonOps_DDRM.invert(A, A_inv);
        if (Double.isInfinite(A_inv.get(0, 0)) || Double.isNaN(A_inv.get(0, 0))) {
            // for singular matrix A, use pseudo-invert instead.
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(A, A_inv);
        }
        CommonOps_DDRM.mult(A_inv, B, phi);
    }

    /**
     * learn parameter Phi with Incremental Computation
     */
    protected void estimateIC() {
        if (-1 != r) {
            DMatrixRMaj subMat = new DMatrixRMaj(dim, dim);
            DMatrixRMaj mat1 = new DMatrixRMaj(dim, dim);
            DMatrixRMaj mat2 = new DMatrixRMaj(dim, dim);
            DMatrixRMaj mat3 = new DMatrixRMaj(dim, dim);
            DMatrixRMaj mat4 = new DMatrixRMaj(dim, dim);
            r += 1;
            for (int i = 1; i <= p; i++) {
                // a_ii
                if (r >= p + 1 - i && r <= tsLength - i) {
                    CommonOps_DDRM.extract(A, (i - 1) * dim, (i - 1) * dim + dim, (i - 1) * dim, (i - 1) * dim + dim, subMat);
                    CommonOps_DDRM.multTransB(ts.getResidualAt(r - 1), ts.getResidualAt(r - 1), mat1);
                    CommonOps_DDRM.multTransB(ts.getLastResidualAt(r - 1), ts.getLastResidualAt(r - 1), mat2);
                    CommonOps_DDRM.subtract(mat1, mat2, mat3);
                    CommonOps_DDRM.add(subMat, mat3, mat1);
                    for (int m = 0; m < dim; m++) {
                        for (int n = 0; n < dim; n++) {
                            A.set((i - 1) * dim + m, (i - 1) * dim + n, mat1.get(m, n));
                        }
                    }

                }
                // a_ij, a_ji, assert i<j
                for (int j = i + 1; j <= p; j++) {
                    if (r >= p + 1 - j && r < p + 1 - i) {
                        CommonOps_DDRM.extract(A, (i - 1) * dim, (i - 1) * dim + dim, (j - 1) * dim, (j - 1) * dim + dim, subMat);
                        CommonOps_DDRM.subtract(ts.getResidualAt(r - 1), ts.getLastResidualAt(r - 1), mat1);
                        CommonOps_DDRM.multTransB(ts.getLastResidualAt(r + j - i - 1), mat1, mat2);
                        CommonOps_DDRM.add(subMat, mat2, mat1);
                        for (int m = 0; m < dim; m++) {
                            for (int n = 0; n < dim; n++) {
                                A.set((i - 1) * dim + m, (j - 1) * dim + n, mat1.get(m, n));
                                A.set((j - 1) * dim + m, (i - 1) * dim + n, mat1.get(n, m));
                            }
                        }
                    } else if (r > tsLength - j && r <= tsLength - i) {
                        CommonOps_DDRM.extract(A, (i - 1) * dim, (i - 1) * dim + dim, (j - 1) * dim, (j - 1) * dim + dim, subMat);
                        CommonOps_DDRM.subtract(ts.getResidualAt(r - 1), ts.getLastResidualAt(r - 1), mat1);
                        CommonOps_DDRM.multTransB(mat1, ts.getLastResidualAt(r - j + i - 1), mat2);
                        CommonOps_DDRM.add(subMat, mat2, mat1);
                        for (int m = 0; m < dim; m++) {
                            for (int n = 0; n < dim; n++) {
                                A.set((i - 1) * dim + m, (j - 1) * dim + n, mat1.get(m, n));
                                A.set((j - 1) * dim + m, (i - 1) * dim + n, mat1.get(n, m));
                            }
                        }
                    } else if (r >= p + 1 - i && r <= tsLength - j) {
                        CommonOps_DDRM.extract(A, (i - 1) * dim, (i - 1) * dim + dim, (j - 1) * dim, (j - 1) * dim + dim, subMat);
                        CommonOps_DDRM.subtract(ts.getResidualAt(r - 1), ts.getLastResidualAt(r - 1), mat1);
                        CommonOps_DDRM.multTransB(ts.getLastResidualAt(r + j - i - 1), mat1, mat2);
                        CommonOps_DDRM.multTransB(mat1, ts.getLastResidualAt(r - j + i - 1), mat3);
                        CommonOps_DDRM.add(mat2, mat3, mat4);
                        CommonOps_DDRM.add(subMat, mat4, mat1);
                        for (int m = 0; m < dim; m++) {
                            for (int n = 0; n < dim; n++) {
                                A.set((i - 1) * dim + m, (j - 1) * dim + n, mat1.get(m, n));
                                A.set((j - 1) * dim + m, (i - 1) * dim + n, mat1.get(n, m));
                            }
                        }
                    } // else r<p+1-j || r>n-i, no change
                }
                // bi
                if (r >= p + 1 - i && r < p + 1) {
                    CommonOps_DDRM.extract(B, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, subMat);
                    CommonOps_DDRM.subtract(ts.getResidualAt(r - 1), ts.getLastResidualAt(r - 1), mat1);
                    CommonOps_DDRM.multTransB(mat1, ts.getLastResidualAt(r + i - 1), mat2);
                    CommonOps_DDRM.add(subMat, mat2, mat1);
                    for (int m = 0; m < dim; m++) {
                        for (int n = 0; n < dim; n++) {
                            B.set((i - 1) * dim + m, n, mat1.get(m, n));
                        }
                    }
                } else if (r > tsLength - i) {
                    CommonOps_DDRM.extract(B, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, subMat);
                    CommonOps_DDRM.subtract(ts.getResidualAt(r - 1), ts.getLastResidualAt(r - 1), mat1);
                    CommonOps_DDRM.multTransB(ts.getLastResidualAt(r - i - 1), mat1, mat2);
                    CommonOps_DDRM.add(subMat, mat2, mat1);
                    for (int m = 0; m < dim; m++) {
                        for (int n = 0; n < dim; n++) {
                            B.set((i - 1) * dim + m, n, mat1.get(m, n));
                        }
                    }
                } else if (r >= p + 1 && r <= tsLength - i) {
                    CommonOps_DDRM.extract(B, (i - 1) * dim, (i - 1) * dim + dim, 0, dim, subMat);
                    CommonOps_DDRM.subtract(ts.getResidualAt(r - 1), ts.getLastResidualAt(r - 1), mat1);
                    CommonOps_DDRM.multTransB(mat1, ts.getLastResidualAt(r + i - 1), mat2);
                    CommonOps_DDRM.multTransB(ts.getLastResidualAt(r - i - 1), mat1, mat3);
                    CommonOps_DDRM.add(mat2, mat3, mat4);
                    CommonOps_DDRM.add(subMat, mat4, mat1);
                    for (int m = 0; m < dim; m++) {
                        for (int n = 0; n < dim; n++) {
                            B.set((i - 1) * dim + m, n, mat1.get(m, n));
                        }
                    }
                } // else r>p+1-i, no change
            }
            ts.setLastRepairAt(r - 1, new DMatrixRMaj(ts.getRepairAt(r - 1)));
            // Restore the r
            r -= 1;
        }
        DMatrixRMaj A_inv = new DMatrixRMaj(dim * p, dim * p);
        CommonOps_DDRM.invert(A, A_inv);
        if (Double.isInfinite(A_inv.get(0, 0)) || Double.isNaN(A_inv.get(0, 0))) {
            // for singular matrix A, use pseudo-invert instead.
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(A, A_inv);
        }
        CommonOps_DDRM.mult(A_inv, B, phi);
    }

    /**
     * Generate a sequence of candidate repair values
     */
    protected void candidate() {
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
                    ts.setCandidateAt(t, new DMatrixRMaj(temp));
                }
            }
        }
    }

    /**
     * Evaluate candidates and perform repairs
     */
    protected void evaluate() {
        double l2norm;
        double minVal = Double.MAX_VALUE;
        int minIndex = -1;
        for (Integer i : ts.candidateList) {
            if (ts.getCandidateAt(i) == null) {
                throw new RuntimeException("[WARNING]null candidate at " + i);
            }
            if (!ts.getLabelAt(i)) {
                l2norm = MyMathUtils.calL2Norm(ts.getCandidateAt(i), ts.getObsAt(i), dim);
                if (l2norm < minVal && l2norm > threshold) {
                    minVal = l2norm;
                    minIndex = i;
                }
            }
        }
        if (minIndex != -1) {
            ts.setRepairAt(minIndex, new DMatrixRMaj(ts.getCandidateAt(minIndex)));
            r = minIndex;
        } else {
            isConverge = true;
        }
        ts.clearCandidates();
    }
}
