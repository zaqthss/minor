package com.MINOR.core;

import com.MINOR.Utils.MyFileUtils;
import com.MINOR.entity.RepairResult;
import org.ejml.data.DMatrixRMaj;

public abstract class BaseCleaningModel {
    protected int dim;
    protected int tsLength;
    protected int p;
    protected double threshold;
    protected String modelName;

    protected DMatrixRMaj Z;
    protected DMatrixRMaj V;
    protected DMatrixRMaj phi;

    BaseCleaningModel(){}

    /**
     * Subclasses need to implement this method to clean time series data and return the repair results.
     *
     * @return A triplet (rmse, timeCost, iterationNum)
     */
    public abstract RepairResult run();

    /**
     * data point with label will be replaced by truth value
     */
    protected abstract void preprocess();

    /**
     * initializing the matrix of OLS, i.e, Z and V
     */
    protected abstract void OLSPre();

    /**
     * learn phi of OLS
     */
    protected abstract void OLS();

    /**
     * calculate the RMSE of the dirty data
     *
     * @return RMSE of the dirty data
     */
    public abstract double getRMSEDirty();

    /**
     * calculate the RMSE of the cleaned data
     *
     * @return RMSE of the cleaned data
     */
    public abstract double getRMSE();

    /**
     * export the time series after cleaning
     *
     * @param fileName file to save the data
     */
    public abstract void exportTimeSeries(String fileName);

    /**
     * export the parameter phi
     *
     * @param fileName file to save the file
     */
    public void exportPhi(String fileName) {
        MyFileUtils.PhiExporter(fileName, this.phi);
    }
}
