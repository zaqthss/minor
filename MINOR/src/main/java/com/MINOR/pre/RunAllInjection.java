package com.MINOR.pre;

/**
 * get all dirty data needed for experiments
 */
public class RunAllInjection {
    public static void main(String[] args) {
        // for error rate, error length, size and error pattern experiment on ILD and ECG dataset
        DataErrInjection.main(args);
        // for gps experiment
        HandleGPS.main(args);
        // for online computing
        ILD48KErrInjection.main(args);
    }
}
