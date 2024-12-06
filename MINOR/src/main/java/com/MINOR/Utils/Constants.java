package com.MINOR.Utils;

public class Constants {
    public final static String dataPathPrefix = "data/";
    public final static String rawPathPrefix = "data/raw/";
    public final static String dirtyPathPrefix = "data/dirty/";
    public final static String appPathPrefix = "data/app/";
    @Deprecated
    public final static String labelPathPrefix = "data/label/";
    public final static String resultPathPrefix = "result/";
    @Deprecated
    public final static String timestampPathPrefix = "data/timestamp/";
    public final static String summaryPathPrefix = "result/summary/";
    public final static String applicationPathPrefix = "D:/MyDocuments/datasets/MINOR-Application/";
    public final static String appCleanPath = "UCRArchive_2018_csv_no_labels/";
    public final static String appDirtyPath = "dirty/";
    public final static String appLabelPath = "label/";

    public final static Long[] seeds = {100L, 200L, 300L, 400L, 500L, 600L, 700L, 800L, 900L, 1000L};
    public final static Long[] seedsApp = {111L, 222L, 777L};

    /**
     * variable setting in experiment
     */
    public final static double[] ers = {0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275};
    public final static int[] els = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
    public final static int[] sizesILD = {5000, 7500, 10000, 12500, 15000, 17500};
    public final static int[] sizesECG = {7000, 8000, 9000, 10000, 11000, 12000,};
    public final static int[] sizesOnline = {5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 48018};
    public final static int[] ps = {1, 2, 3, 4, 5, 6};
    public final static double[] lrs = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
    public final static double[] thresholds = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};
    public final static int[] maxIterations = {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000};
    //
    public final static double er0 = 0.05;
    public final static int el0 = 50;
    public final static int p0 = 1;
    public final static double lr0 = 0.2;
    public final static double thresholdGps0 = 0.05;
    public final static double thresholdILD0 = 0.05;
    public final static double thresholdECG0 = 0.02;
    public final static int maxIteration0 = 100_000;
}

