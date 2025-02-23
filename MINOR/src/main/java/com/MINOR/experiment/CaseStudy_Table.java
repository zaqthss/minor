package com.MINOR.experiment;

import com.MINOR.Utils.Constants;
import com.MINOR.Utils.MyFileUtils;
import com.MINOR.Utils.MyMathUtils;
import com.MINOR.core.*;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimeSeries;

import java.util.Arrays;
import java.util.List;

import static java.lang.System.exit;

/**
 * result of TABLE IV
 */
public class CaseStudy_Table {
    public static void main(String[] args) {
        // experiment setting
        int dim = 2;
        int p0 = Constants.p0;
        double lr0 = Constants.lr0;
        double[] thresholds = {0.05};
        int maxIterations = Constants.maxIteration0;
        int caseCnt = thresholds.length;
        int repeatTime = Constants.seeds.length;
        double S = 1.4325;      // speed constraint of MTCSC-C
        int T = 100;            // window size of MTCSC-C
        double SMax = 1.4325;   // speed constraint of MTCSC-A
        int K = 3;              // order of Markov-chain for Akane Heuristic Auto
        double tolerance = 0.000001; // eliminate precision issues
        // results to export
        String[] modelNames = {"Dirty", "MINOR-B", "MINOR-B-w/o CIC", "MINOR-B-w/o IC", "MINOR-B-Uni", "MINOR-B-w/o SC", "MINOR-B-w/o iter."
                , "MINOR-U", "MINOR-U-w/o CIC", "MINOR-U-w/o IC", "MINOR-U-w/o SC", "MINOR-U-w/o iter."
                , "VARX", "IMR", "MTCSC-C", "MTCSC-A", "Akane"};
        int modelCnt = modelNames.length;
        double[][] totalRMSE = new double[caseCnt][modelCnt];
        double[][] totalTime = new double[caseCnt][modelCnt];
        int[][] totalIterNum = new int[caseCnt][modelCnt];
        // local variables in loop
        TimeSeries ts;
        List<UniTimeSeries> unis;
        RepairResult result;
        String dataSetName;
        String dataFPath;
        int modelID;
        double[] RMSEOfDim = new double[dim];
        double[] timeOfDim = new double[dim];
        int[] iterNumOfDim = new int[dim];
        // determine whether each model will be used
        boolean[] flags = new boolean[modelCnt];
        Arrays.fill(flags, true);

        for (int i = 0; i < caseCnt; i++) {
            for (int seed = 0; seed < repeatTime; seed++) {
                ts = new TimeSeries();
                dataSetName = "GPS/GPS_lr_" + lr0 + '_' + seed;
                dataFPath = Constants.dirtyPathPrefix + dataSetName + ".csv";
                System.out.println("--------------------");
                System.out.println("[INFO]testing " + dataSetName);
                double threshold = thresholds[i];
                // Dirty
                modelID = 0;
                MyFileUtils.csv2MTS(ts, dim, dataFPath);
                totalRMSE[i][modelID] += MyMathUtils.calculateRMSE_dirty(ts);
                // elimate the starting lag
                MyFileUtils.csv2MTS(ts, dim, dataFPath);
                MINOR_B minorB2 = new MINOR_B(ts, p0, threshold, maxIterations);
                result = minorB2.run();
                // MINOR-B
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_B minorB = new MINOR_B(ts, p0, threshold, maxIterations);
                    result = minorB.run();
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-B-w/o CIC
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_B minorB = new MINOR_B(ts, p0, threshold, maxIterations);
                    RepairResult result2 = minorB.runWithoutCandiIC();
                    if (result2.getRmse() - totalRMSE[i][modelID - 1] > tolerance) {
                        System.out.println("[WARNING]Exceeded the precision limit!");
                        exit(1);
                    }
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result2.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-B-w/o IC
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_B minorB = new MINOR_B(ts, p0, threshold, maxIterations);
                    RepairResult result2 = minorB.runSD();
                    if (result2.getRmse() - totalRMSE[i][modelID - 1] > tolerance) {
                        System.out.println("[WARNING]Exceeded the precision limit!");
                        exit(1);
                    }
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result2.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-B-Uni
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    unis = ts.toUniTimeSeries();
                    for (int m = 0; m < dim; m++) {
                        MINOR_BUni minor_bUni = new MINOR_BUni(unis.get(m), p0, threshold, maxIterations);
                        result = minor_bUni.run();
                        RMSEOfDim[m] = result.getRmse();
                        timeOfDim[m] = result.getTimeCost();
                        iterNumOfDim[m] = result.getIterationNum();
                    }
                    double rmseTmp = 0;
                    double timeTmp = 0;
                    int iterTmp = Arrays.stream(iterNumOfDim).max().orElseThrow();
                    for (int m = 0; m < dim; m++) {
                        rmseTmp += RMSEOfDim[m] * RMSEOfDim[m];
                        timeTmp += timeOfDim[m];
                    }
                    rmseTmp = Math.sqrt(rmseTmp);
                    totalRMSE[i][modelID] += rmseTmp;
                    totalTime[i][modelID] += timeTmp;
                    totalIterNum[i][modelID] += iterTmp;
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + new RepairResult(rmseTmp, timeTmp, iterTmp) + '\n');
                }
                // MINOR-B-w/o SC
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_B_no_sc minorB = new MINOR_B_no_sc(ts, p0, threshold, maxIterations);
                    result = minorB.run();
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-B-w/o iter.
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_B_no_iter minorB = new MINOR_B_no_iter(ts, p0, threshold);
                    result = minorB.run();
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-U
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_U minorU = new MINOR_U(ts, p0, threshold, maxIterations);
                    result = minorU.run();
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-U-w/o CIC
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_U minorU = new MINOR_U(ts, p0, threshold, maxIterations);
                    RepairResult result2 = minorU.runWithoutCandiIC();
                    if (result2.getRmse() - totalRMSE[i][modelID - 1] > tolerance) {
                        System.out.println("[WARNING]Exceeded the precision limit!");
                        exit(1);
                    }
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result2.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-U-w/o IC
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_U minorU = new MINOR_U(ts, p0, threshold, maxIterations);
                    RepairResult result2 = minorU.runSD();
                    if (result2.getRmse() - totalRMSE[i][modelID - 1] > tolerance) {
                        System.out.println("[WARNING]Exceeded the precision limit!");
                        exit(1);
                    }
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result2.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-U-w/o SC
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_base minorU = new MINOR_base(ts, p0, threshold, maxIterations);
                    result = minorU.run();
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // MINOR-U-w/o iter.
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_U_no_iter minorU = new MINOR_U_no_iter(ts, p0, threshold);
                    result = minorU.run();
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // VARX
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    VARX varx = new VARX(ts, p0, threshold);
                    result = varx.run();
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // IMR
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    unis = ts.toUniTimeSeries();
                    for (int m = 0; m < dim; m++) {
                        IMR imr = new IMR(unis.get(m), p0, threshold, maxIterations);
                        result = imr.run();
                        RMSEOfDim[m] = result.getRmse();
                        timeOfDim[m] = result.getTimeCost();
                        iterNumOfDim[m] = result.getIterationNum();
                    }
                    double rmseTmp = 0;
                    double timeTmp = 0;
                    int iterTmp = Arrays.stream(iterNumOfDim).max().orElseThrow();
                    for (int m = 0; m < dim; m++) {
                        rmseTmp += RMSEOfDim[m] * RMSEOfDim[m];
                        timeTmp += timeOfDim[m];
                    }
                    rmseTmp = Math.sqrt(rmseTmp);
                    totalRMSE[i][modelID] += rmseTmp;
                    totalTime[i][modelID] += timeTmp;
                    totalIterNum[i][modelID] += iterTmp;
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + new RepairResult(rmseTmp, timeTmp, iterTmp) + '\n');
                }
                // MTCSC-C
                modelID++;
                if (flags[modelID]) {
                    if (i > 0) {
                        totalRMSE[i][modelID] += totalRMSE[0][modelID];
                        totalTime[i][modelID] += totalTime[0][modelID];
                        totalIterNum[i][modelID] += totalIterNum[0][modelID];
                    } else {
                        System.out.println("[INFO]running " + modelNames[modelID]);
                        MyFileUtils.csv2MTS(ts, dim, dataFPath);
                        MTCSC mtcsc_c = new MTCSC(ts, S, T, dim);
                        result = mtcsc_c.mainScreen();
                        totalRMSE[i][modelID] += result.getRmse();
                        totalTime[i][modelID] += result.getTimeCost();
                        totalIterNum[i][modelID] += result.getIterationNum();
                        System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                    }
                }
                // MTCSC-A
                modelID++;
                if (flags[modelID]) {
                    if (i > 0) {
                        totalRMSE[i][modelID] += totalRMSE[0][modelID];
                        totalTime[i][modelID] += totalTime[0][modelID];
                        totalIterNum[i][modelID] += totalIterNum[0][modelID];
                    } else {
                        System.out.println("[INFO]running " + modelNames[modelID]);
                        MyFileUtils.csv2MTS(ts, dim, dataFPath);
                        MTCSC_A mtcsc_a = new MTCSC_A(ts, SMax, T, 0.05, 0.75, 125, 0.75, 5);
                        result = mtcsc_a.mainScreen();
                        totalRMSE[i][modelID] += result.getRmse();
                        totalTime[i][modelID] += result.getTimeCost();
                        totalIterNum[i][modelID] += result.getIterationNum();
                        System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                    }
                }
                // Akane
                modelID++;
                if (flags[modelID]) {
                    if (i > 0) {
                        totalRMSE[i][modelID] += totalRMSE[0][modelID];
                        totalTime[i][modelID] += totalTime[0][modelID];
                        totalIterNum[i][modelID] += totalIterNum[0][modelID];
                    } else {
                        System.out.println("[INFO]running " + modelNames[modelID]);
                        MyFileUtils.csv2MTS(ts, dim, dataFPath);
                        unis = ts.toUniTimeSeries();
                        for (int m = 0; m < dim; m++) {
                            AkaneHeuristic akaneheuristicauto = new AkaneHeuristic(unis.get(m), K);
                            result = akaneheuristicauto.mainAkaneHeuristicAuto();
                            RMSEOfDim[m] = result.getRmse();
                            timeOfDim[m] = result.getTimeCost();
                            iterNumOfDim[m] = result.getIterationNum();
                        }
                        double rmseTmp = 0;
                        double timeTmp = 0;
                        int iterTmp = 1;
                        for (int m = 0; m < dim; m++) {
                            rmseTmp += RMSEOfDim[m] * RMSEOfDim[m];
                            timeTmp += timeOfDim[m];
                        }
                        rmseTmp = Math.sqrt(rmseTmp);
                        totalRMSE[i][modelID] += rmseTmp;
                        totalTime[i][modelID] += timeTmp;
                        totalIterNum[i][modelID] += iterTmp;
                        System.out.println("[RESULT]" + modelNames[modelID] + ":" + new RepairResult(rmseTmp, timeTmp, iterTmp) + '\n');
                    }
                }
            }
            for (int j = 0; j < modelCnt; j++) {
                totalRMSE[i][j] /= repeatTime;
                totalTime[i][j] /= repeatTime;
                totalIterNum[i][j] /= repeatTime;
            }
        }
        String basePath = Constants.resultPathPrefix + "CaseStudy/table/";
        String fileName = basePath + "RMS.csv";
        MyFileUtils.exportResult(fileName, modelNames, thresholds, totalRMSE);
        fileName = basePath + "TIME.csv";
        MyFileUtils.exportResult(fileName, modelNames, thresholds, totalTime);
        fileName = basePath + "iterationNum.csv";
        MyFileUtils.exportResult(fileName, modelNames, thresholds, totalIterNum);
    }
}
