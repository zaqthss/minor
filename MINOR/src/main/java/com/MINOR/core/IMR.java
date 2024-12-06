package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.UniTimeSeries;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.ops.SimpleOperations_DDRM;

public class IMR extends UTSCModel {
    protected long iterationLimit;
    protected boolean isConverge;

    protected DMatrixRMaj A;
    protected DMatrixRMaj B;
    protected int r;

    public IMR(UniTimeSeries ts, int p, double threshold, long iterationLimit) {
        this.modelName = "IMR";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;
        this.iterationLimit = iterationLimit;
        r = -1;
        this.isConverge = false;
        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, p);
        V = new DMatrixRMaj(tsLength - p, 1);
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        this.phi = new DMatrixRMaj(p, 1);
    }

    public IMR() {
    }

    @Override
    public RepairResult run(){
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
     * used for application experiment
     *
     * @return cleaned data
     */
    public UniTimeSeries runApp() {
        preprocess();
        MatrixPruning();
        int k = 1;
        this.isConverge = false;
        while (true) {
            estimateIC();
            candidate();
            evaluate();
            if (isConverge || k >= iterationLimit) {
                break;
            }
            k++;
        }
        return ts;
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
        for (int i = 0; i < tsLength - p; i++) {
            for (int j = 0; j < p; j++) {
                Z.set(i, j, ts.getResidualAt(p - 1 + i - j));
            }
        }
        for (int i = 0; i < tsLength - p; i++) {
            V.set(i, 0, ts.getResidualAt(p + i));
        }
    }

    /**
     * remove rows with all zeros of Z and the corresponding rows of V, and generate matrix A and B for OLS
     */
    protected void MatrixPruning() {
        OLSPre();
        DMatrixRMaj[] wrapper = {Z, V};
        MyMathUtils.MatrixPruning(wrapper);
        Z = wrapper[0];
        V = wrapper[1];
        A.reshape(Z.numCols, Z.numCols);
        B.reshape(Z.numCols, V.numCols);
        CommonOps_DDRM.multTransA(Z, Z, A);         // A_0 = Z'*Z
        CommonOps_DDRM.multTransA(Z, V, B);         // B_0 = Z'*V
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
        DMatrixRMaj A_inv = new DMatrixRMaj(p, p);
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
            double subMat, mat1, mat2, mat3, mat4;
            r += 1;                   // Align the range of values for r with the r in the paper.
            for (int i = 1; i <= p; i++) {
                // a_ii
                if (r >= p + 1 - i && r <= tsLength - i) {
                    subMat = A.get(i - 1, i - 1);
                    mat1 = ts.getResidualAt(r - 1) * ts.getResidualAt(r - 1);
                    mat2 = ts.getLastResidualAt(r - 1) * ts.getLastResidualAt(r - 1);
                    mat3 = mat1 - mat2;
                    mat1 = subMat + mat3;
                    A.set((i - 1), (i - 1), mat1);
                }
                // a_ij, a_ji, assert i<j
                for (int j = i + 1; j <= p; j++) {
                    if (r >= p + 1 - j && r < p + 1 - i) {
                        subMat = A.get(i - 1, j - 1);
                        mat1 = ts.getResidualAt(r - 1) - ts.getLastResidualAt(r - 1);
                        mat2 = ts.getLastResidualAt(r + j - i - 1) * mat1;
                        mat1 = subMat + mat2;
                        A.set(j - 1, i - 1, mat1);
                        A.set(i - 1, j - 1, mat1);
                    } else if (r > tsLength - j && r <= tsLength - i) {
                        subMat = A.get(i - 1, j - 1);
                        mat1 = ts.getResidualAt(r - 1) - ts.getLastResidualAt(r - 1);
                        mat2 = ts.getLastResidualAt(r - j + i - 1) * mat1;
                        mat1 = subMat + mat2;
                        A.set(j - 1, i - 1, mat1);
                        A.set(i - 1, j - 1, mat1);
                    } else if (r >= p + 1 - i && r <= tsLength - j) {
                        subMat = A.get(i - 1, j - 1);
                        mat1 = ts.getResidualAt(r - 1) - ts.getLastResidualAt(r - 1);
                        mat2 = ts.getLastResidualAt(r + j - i - 1) * mat1;
                        mat3 = ts.getLastResidualAt(r - j + 1 - 1) * mat1;
                        mat4 = mat2 + mat3;
                        mat1 = subMat + mat4;
                        A.set(j - 1, i - 1, mat1);
                        A.set(i - 1, j - 1, mat1);
                    } // else r<p+1-j || r>n-i, no change
                }
                // bi
                if (r >= p + 1 - i && r < p + 1) {
                    subMat = B.get(i - 1, 0);
                    mat1 = ts.getResidualAt(r - 1) - ts.getLastResidualAt(r - 1);
                    mat2 = mat1 * ts.getLastResidualAt(r + i - 1);
                    mat1 = subMat + mat2;
                    B.set(i - 1, 0, mat1);
                } else if (r > tsLength - i) {
                    subMat = B.get(i - 1, 0);
                    mat1 = ts.getResidualAt(r - 1) - ts.getLastResidualAt(r - 1);
                    mat2 = mat1 * ts.getLastResidualAt(r - i - 1);
                    mat1 = subMat + mat2;
                    B.set(i - 1, 0, mat1);
                } else if (r >= p + 1 && r <= tsLength - i) {
                    subMat = B.get(i - 1, 0);
                    mat1 = ts.getResidualAt(r - 1) - ts.getLastResidualAt(r - 1);
                    mat2 = mat1 * ts.getLastResidualAt(r + i - 1);
                    mat3 = ts.getLastResidualAt(r - i - 1) * mat1;
                    mat4 = mat2 + mat3;
                    mat1 = subMat + mat4;
                    B.set(i - 1, 0, mat1);
                } // else r>p+1-i, no change
            }
            ts.setLastRepairAt(r - 1, ts.getRepairAt(r - 1));
            // Restore the r
            r -= 1;
        }
        DMatrixRMaj A_inv = new DMatrixRMaj(p, p);
        CommonOps_DDRM.invert(A, A_inv);
        if (Double.isInfinite(A_inv.get(0, 0)) || Double.isNaN(A_inv.get(0, 0))) {
            SimpleOperations_DDRM simpleOperations_ddrm = new SimpleOperations_DDRM();
            simpleOperations_ddrm.pseudoInverse(A, A_inv);
        }
        CommonOps_DDRM.mult(A_inv, B, phi);
    }

    /**
     * Generate a sequence of candidate repair values
     */
    private void candidate() {
        double temp;
        Boolean lb;
        for (int t = p; t < tsLength; t++) {
            lb = ts.getLabelAt(t);
            if (!lb) {
                temp = 0;
                for (int i = 1; i <= p; i++) {
                    temp += phi.get(i - 1, 0) * ts.getResidualAt(t - i);
                }
                temp += ts.getObsAt(t);
                double l2norm = Math.abs(temp - ts.getRepairAt(t));
                if (l2norm > this.threshold) {
                    ts.setCandidateAt(t, temp);
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
        boolean lb;
        for (Integer i : ts.candidateList) {
            if (ts.getCandidateAt(i) == Double.MAX_VALUE) {
                throw new RuntimeException("[WARNING]null candidate at " + i);
            }
            lb = ts.getLabelAt(i);
            if (!lb) {
                l2norm = Math.abs(ts.getCandidateAt(i) - ts.getObsAt(i));
                if (l2norm < minVal && l2norm > threshold) {
                    minVal = l2norm;
                    minIndex = i;
                }
            }
        }
        if (minIndex != -1) {
            ts.setRepairAt(minIndex, ts.getCandidateAt(minIndex));
            r = minIndex;
        } else {
            isConverge = true;
        }
        ts.clearCandidates();
    }
}
