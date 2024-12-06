package com.MINOR.experiment;

import com.MINOR.Utils.Constants;
import com.MINOR.Utils.MyFileUtils;
import com.MINOR.core.*;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimeSeries;

import java.util.Arrays;
import java.util.List;

/**
 * Varying maximum iteration number on GPS
 */
public class MaxIteration_GPS {
    public static void main(String[] args) {
        // experiment setting
        int dim = 2;
        String variableName = "maxiteration";
        int p0 = Constants.p0;
        double lr0 = Constants.lr0;
//        double threshold0 = Constants.thresholdGps0;
        double threshold0 = 0.02;
        int[] maxIterations = Constants.maxIterations;
        int caseCnt = maxIterations.length;
        int repeatTime = Constants.seeds.length;
        double S = 1.4325;      // speed constraint of MTCSC-C
        int T = 100;            // window size of MTCSC-C
        double SMax = 1.4325;   // speed constraint of MTCSC-A
        int K = 3;              // order of Markov-chain for Akane Heuristic Auto
        // results to export
        String[] modelNames = {"MINOR-B", "MINOR-U", "VARX", "IMR", "MTCSC-C", "Akane", "MTCSC-A"};
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
                int maxIteration = maxIterations[i];
                // MINOR-B
                modelID = 0;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_B minorB = new MINOR_B(ts, p0, threshold0, maxIteration);
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
                    MINOR_U minorU = new MINOR_U(ts, p0, threshold0, maxIteration);
                    result = minorU.run();
                    totalRMSE[i][modelID] += result.getRmse();
                    totalTime[i][modelID] += result.getTimeCost();
                    totalIterNum[i][modelID] += result.getIterationNum();
                    System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                }
                // VARX
                modelID++;
                if (flags[modelID]) {
                    if (i > 0) {
                        totalRMSE[i][modelID] += totalRMSE[0][modelID];
                        totalTime[i][modelID] += totalTime[0][modelID];
                        totalIterNum[i][modelID] += totalIterNum[0][modelID];
                    } else {
                        System.out.println("[INFO]running " + modelNames[modelID]);
                        MyFileUtils.csv2MTS(ts, dim, dataFPath);
                        VARX varx = new VARX(ts, p0, threshold0);
                        result = varx.run();
                        totalRMSE[i][modelID] += result.getRmse();
                        totalTime[i][modelID] += result.getTimeCost();
                        totalIterNum[i][modelID] += result.getIterationNum();
                        System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');
                    }
                }
                // IMR
                modelID++;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    unis = ts.toUniTimeSeries();
                    for (int m = 0; m < dim; m++) {
                        IMR imr = new IMR(unis.get(m), p0, threshold0, maxIteration);
                        result = imr.run();
                        RMSEOfDim[m] = result.getRmse();
                        timeOfDim[m] = result.getTimeCost();
                        iterNumOfDim[m] = result.getIterationNum();
                    }
                    double rmseTmp = 0;
                    double timeTmp = 0;
                    int iterTmp = 0;
                    for (int m = 0; m < dim; m++) {
                        rmseTmp += RMSEOfDim[m] * RMSEOfDim[m];
                        timeTmp += timeOfDim[m];
                        iterTmp += iterNumOfDim[m];
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
            }
            for (int j = 0; j < modelCnt; j++) {
                totalRMSE[i][j] /= repeatTime;
                totalTime[i][j] /= repeatTime;
                totalIterNum[i][j] /= repeatTime;
            }
        }
        String basePath = Constants.resultPathPrefix + "GPS/summary/" + variableName + '/';
        String fileName = basePath + "RMS.csv";
        MyFileUtils.exportResult(fileName, modelNames, maxIterations, totalRMSE);
        fileName = basePath + "TIME.csv";
        MyFileUtils.exportResult(fileName, modelNames, maxIterations, totalTime);
        fileName = basePath + "iterationNum.csv";
        MyFileUtils.exportResult(fileName, modelNames, maxIterations, totalIterNum);
    }
}
