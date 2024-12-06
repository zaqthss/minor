package com.MINOR.experiment;

import com.MINOR.Utils.Constants;
import com.MINOR.Utils.MyFileUtils;
import com.MINOR.core.*;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimeSeries;
import com.MINOR.enums.DataSetType;
import com.MINOR.enums.ErrorInjectorType;

import java.util.Arrays;
import java.util.List;

public class OnlineComputing {
    public static void main(String[] args) {
//        run(DataSetType.ILD, ErrorInjectorType.SHIFT);

        // experiment setting
        int dim = 2;
        String variableName = "size";
        DataSetType dt = DataSetType.ILD;
        ErrorInjectorType et = ErrorInjectorType.SHIFT;
        int p0 = Constants.p0;
        double threshold0 = Constants.thresholdILD0;
        int maxIterations = Constants.maxIteration0;
        int[] sizes = Constants.sizesOnline;
        int caseCnt = sizes.length;
        int repeatTime = Constants.seeds.length * 10;
        double S = 0.06825;      // speed constraint of MTCSC-C
        int T = 100;            // window size of MTCSC-C
        int K = 3;              // order of Markov-chain for Akane Heuristic Auto
        int labelLen = 1000;    // for IMR stream
        // results to export
        String[] modelNames = {"MINOR-O", "VARX", "IMR", "MTCSC"};
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
                dataSetName = "ILD48k_ONLINE_" + variableName + "_" + sizes[i] + '_' + seed % 10 + ".csv";
                dataFPath = Constants.dirtyPathPrefix + "ILD48k/" + dataSetName;
                System.out.println("--------------------");
                System.out.println("[INFO]testing " + dataSetName);
                // MINOR-O
                modelID=0;
                if (flags[modelID]) {
                    System.out.println("[INFO]running " + modelNames[modelID]);
                    MyFileUtils.csv2MTS(ts, dim, dataFPath);
                    MINOR_O minorO = new MINOR_O(ts);
                    result = minorO.run();
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
                    VARXStream varxs = new VARXStream(ts, threshold0);
                    result = varxs.run();
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
                        IMRStream imrs = new IMRStream(unis.get(m), labelLen);
                        result = imrs.run();
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
                // MTCSC-C
                modelID++;
                if (flags[modelID]) {
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
            for (int j = 0; j < modelCnt; j++) {
                totalRMSE[i][j] /= repeatTime;
                totalTime[i][j] /= repeatTime;
                totalIterNum[i][j] /= repeatTime;
            }
        }
        double[][] totalTimeDiff = new double[caseCnt - 1][modelCnt];
        double[][] totalThroughout = new double[caseCnt - 1][modelCnt];
        for (int j = 0; j < modelCnt; j++) {
            for (int i = 0; i < caseCnt - 1; i++) {
                totalTimeDiff[i][j] = totalTime[i + 1][j] - totalTime[i][j];
                if (totalTimeDiff[i][j] == 0) {
                    totalThroughout[i][j] = 0;
                } else {
                    totalThroughout[i][j] = (((double) (sizes[i + 1] - sizes[i])) / totalTimeDiff[i][j]) * 1000 / 1e7;
                }
            }
        }
        double[][] totalRMSE2 = new double[caseCnt - 3][modelCnt];
        double[][] totalThroughout2 = new double[caseCnt - 3][modelCnt];
        double[][] totalIterNum2 = new double[caseCnt - 3][modelCnt];
        double[][] totalTimeDiff2 = new double[caseCnt - 3][modelCnt];
        int[] sizes2 = new int[caseCnt - 3];
        for (int i = 0; i < caseCnt - 3; i++) {
            sizes2[i] = sizes[i + 3];
            for (int j = 0; j < modelCnt; j++) {
                totalRMSE2[i][j] = totalRMSE[i + 3][j];
                totalThroughout2[i][j] = totalThroughout[i + 2][j];
                totalTimeDiff2[i][j] = totalTimeDiff[i + 2][j];
                totalIterNum2[i][j] = totalIterNum[i + 3][j];
            }
        }
        String basePath = Constants.resultPathPrefix + "online/";
        String fileName = basePath + "RMS.csv";
        MyFileUtils.exportResult(fileName, modelNames, sizes2, totalRMSE2);
        fileName = basePath + "throughout.csv";
        MyFileUtils.exportResult(fileName, modelNames, sizes2, totalThroughout2);
//        fileName = basePath + "iterationNum.csv";
//        MyFileUtils.exportResult(fileName, modelNames, sizes2, totalIterNum2);
//        fileName = basePath + "TIME.csv";
//        MyFileUtils.exportResult(fileName, modelNames, sizes, totalTime);
//        fileName = basePath + "TIME_diff.csv";
//        MyFileUtils.exportResult(fileName, modelNames, sizes2, totalTimeDiff2);
    }
}
