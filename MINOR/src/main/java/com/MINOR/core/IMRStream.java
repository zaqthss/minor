package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.UniTimePoint;
import com.MINOR.entity.UniTimeSeries;
import org.ejml.data.DMatrixRMaj;

public class IMRStream extends IMR{
    private final int labelLen;

    public IMRStream(UniTimeSeries ts, int labelLen) {
        this.modelName = "IMRStream";

        this.ts = ts;
        this.p = 1;
        this.threshold = 0;
        this.iterationLimit = 1;
        this.labelLen = labelLen;
        r = -1;
        this.isConverge = false;
        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, p);
        V = new DMatrixRMaj(tsLength - p, 1);
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        this.phi = new DMatrixRMaj(p, 1);
    }

    @Override
    public RepairResult run() {
        if (p != 1) {
            System.out.println("Can not process IMRstream p = " + p);
        }
        long startTime = System.nanoTime();
        // form z
        double[] zs = new double[tsLength];
        UniTimePoint tp;
        for (int i = 0; i < tsLength; ++i) {
            tp = ts.getP(i);
            zs[i] = tp.getLabel() ? tp.getValTruth() - tp.getValObs() : 0;
        }
        double alpha = 0, beta = 0;
        for (int i = 1; i < labelLen - 1; ++i) {
            alpha += zs[i] * zs[i];
            beta += zs[i] * zs[i - 1];
        }
        alpha += zs[0] * zs[0];
        beta += zs[labelLen - 1] * zs[labelLen - 2];
        double phi;
        // learn the original phi
        phi = beta / alpha;
        double estimate;
        boolean isUpdate;
        double val, fillVal;
        // stream input == become not 0 one by one
        for (int readIndex = labelLen; readIndex < tsLength; ++readIndex) {
            tp = ts.getP(readIndex);
            val = zs[readIndex]; // original y-x
            estimate = 0;
            if (tp.getLabel()) {
                isUpdate = false;
            } else {
                estimate = phi * zs[readIndex - 1];
                isUpdate = true;
            }
            // update
            if (isUpdate) {
                fillVal = estimate;
            } else {
                fillVal = val;
            }
            zs[readIndex] = fillVal;
            alpha += zs[readIndex - 1] * zs[readIndex - 1];
            beta += zs[readIndex] * zs[readIndex - 1];
            phi = beta / alpha;
        }
        // begin repairing
        double modify;
        for (int i = 0; i < tsLength; ++i) {
            tp = ts.getP(i);
            if (tp.getLabel()) {
                modify = tp.getValTruth();
            } else {
                modify = tp.getValObs() + zs[i];
            }

            tp.setValRepaired(modify);
        }
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        return new RepairResult(MyMathUtils.calculateRMSE(ts), duration, 1);
    }
}
