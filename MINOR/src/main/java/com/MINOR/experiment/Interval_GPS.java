package com.MINOR.experiment;

import com.MINOR.Utils.Constants;
import com.MINOR.Utils.MyFileUtils;
import com.MINOR.Utils.MyMathUtils;
import com.MINOR.core.*;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimePoint;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimeSeries;

import java.util.*;
import java.util.stream.Collectors;

/**
 * result of the experiment on the interval between labeled and repaired data points
 */
public class Interval_GPS {
    public static void main(String[] args) {
        // experiment setting
        int dim = 2;
        int p0 = Constants.p0;
        double lr0 = Constants.lr0;
        double threshold_0 = 0.09;
        int maxIterations = Constants.maxIteration0;

        int repeatTime = Constants.seeds.length;
        // results to export
        String[] modelNames = {"MINOR-B", "MINOR-U"};
        int modelCnt = modelNames.length;
        // local variables in loop
        TimeSeries ts;
        List<UniTimeSeries> unis;
        RepairResult result;
        String dataSetName;
        String dataFPath;
        int modelID;

        int[] intervalCases = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        int intervalLimit = intervalCases.length;

        List<ErrorInfo> errors_B = new ArrayList<>();
        List<ErrorInfo> errors_U = new ArrayList<>();
        Map<Long, Double> map_U = new HashMap<>();
        Map<Long, Double> map_B = new HashMap<>();
        // 记录interval为特定值的数据点的数量
        int[][] totalNum = new int[intervalLimit][modelCnt];
        double[][] totalRMSE = new double[intervalLimit][modelCnt];
        // 用于避免有些种子里，interval为特定值的数据点数量为0造成的rmse为0带来的最终RMSE的不准确
        int[][] nonZeroCases = new int[intervalLimit][modelCnt];

        for (int seed = 0; seed < repeatTime; seed++) {
            ts = new TimeSeries();
            dataSetName = "GPS/GPS_lr_" + lr0 + '_' + seed;
            dataFPath = Constants.dirtyPathPrefix + dataSetName + ".csv";
            System.out.println("--------------------");
            System.out.println("[INFO]testing " + dataSetName);
            // MINOR-B
            modelID = 0;
            System.out.println("[INFO]running " + modelNames[modelID]);
            MyFileUtils.csv2MTS(ts, dim, dataFPath);
            MINOR_B minorB = new MINOR_B(ts, p0, threshold_0, maxIterations);
            result = minorB.run();
            System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');

            int length = ts.size();
            int lastLabel;
            int nextLabel = -1;
            long interval;
            double err_obs, err_repaired, ratio;
            errors_B.clear();
            map_B.clear();
            lastLabel = -1;
            for (int t = 0; t < length; t++) {
                TimePoint tp = ts.getP(t);
                if (tp.getLabel()) {
                    lastLabel = t;
                } else {
                    if (lastLabel == -1) {
                        continue;
                    }
                    if (tp.isDirty()) {
                        for (int k = t + 1; k < length; k++) {
                            if (ts.getLabelAt(k)) {
                                nextLabel = k;
                                break;
                            }
                            if (k == length - 1) {
                                nextLabel = -1;
                                break;
                            }
                        }
                        if (nextLabel == -1) {
                            interval = tp.getTimestamp() - ts.getTimestampAt(lastLabel);
                        } else {
                            long interval_left = tp.getTimestamp() - ts.getTimestampAt(lastLabel);
                            long interval_right = ts.getTimestampAt(nextLabel) - tp.getTimestamp();
                            interval = Math.min(interval_left, interval_right);
                        }
                        if (interval > intervalLimit) {
                            continue;
                        }
                        interval--;
                        totalNum[(int) interval][modelID]++;
                        err_obs = MyMathUtils.calL2Norm(tp.getValObs(), tp.getValTruth(), dim);
                        err_repaired = MyMathUtils.calL2Norm(tp.getValRepaired(), tp.getValTruth(), dim);
                        ErrorInfo errInfo = new ErrorInfo(interval, err_obs, err_repaired);
                        errors_B.add(errInfo);
                    }
                }
            }

            map_B = calculateRatioMap(errors_B);
            for (int i = 0; i < intervalLimit; i++) {
                if (map_B.get((long) i) == null) {
                    ratio = 0;
                } else {
                    ratio = map_B.get((long) i);
                }
                totalRMSE[i][modelID] += ratio;
            }

            // MINOR-U
            modelID++;
            System.out.println("[INFO]running " + modelNames[modelID]);
            MyFileUtils.csv2MTS(ts, dim, dataFPath);
            MINOR_U minorU = new MINOR_U(ts, p0, threshold_0, maxIterations);
            result = minorU.run();
            System.out.println("[RESULT]" + modelNames[modelID] + ":" + result + '\n');

            length = ts.size();
            errors_U.clear();
            map_U.clear();
            lastLabel = -1;
            for (int t = 0; t < length; t++) {
                TimePoint tp = ts.getP(t);
                if (tp.getLabel()) {
                    lastLabel = t;
                } else {
                    if (lastLabel == -1) {
                        continue;
                    }
                    if (tp.isDirty()) {
                        interval = tp.getTimestamp() - ts.getTimestampAt(lastLabel);
                        if (interval > intervalLimit) {
                            continue;
                        }
                        interval--;
                        totalNum[(int) interval][modelID]++;
                        err_obs = MyMathUtils.calL2Norm(tp.getValObs(), tp.getValTruth(), dim);
                        err_repaired = MyMathUtils.calL2Norm(tp.getValRepaired(), tp.getValTruth(), dim);
                        ErrorInfo errInfo = new ErrorInfo(interval, err_obs, err_repaired);
                        errors_U.add(errInfo);
                    }
                }
            }
            map_U = calculateRatioMap(errors_U);
            for (int i = 0; i < intervalLimit; i++) {
                if (map_U.get((long) i) == null) {
                    ratio = 0;
                } else {
                    ratio = map_U.get((long) i);
                    nonZeroCases[i][0] ++;
                    nonZeroCases[i][1] ++;
                }
                totalRMSE[i][modelID] += ratio;
            }
        }

        for (int i = 0; i < intervalLimit; i++) {
            for (int j = 0; j < modelCnt; j++) {
                totalNum[i][j] /= nonZeroCases[i][j];
                totalRMSE[i][j] /= nonZeroCases[i][j];
            }
        }

        String basePath = Constants.resultPathPrefix + "GPS/summary/interval/";
        String fileName = basePath + "RMSE_ratio.csv";
        MyFileUtils.exportResult(fileName, modelNames, intervalCases, totalRMSE);
        fileName = basePath + "NUM.csv";
        MyFileUtils.exportResult(fileName, modelNames, intervalCases, totalNum);
    }

    private static Map<Long, Double> calculateRatioMap(List<ErrorInfo> ls) {
        return ls.stream()
                .collect(Collectors.groupingBy(e -> e.interval))
                .entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        entry -> {
                            List<ErrorInfo> group = entry.getValue();
                            int count = group.size();

                            // 计算err_obs的RMSE
                            double rmseObs = Math.sqrt(
                                    group.stream()
                                            .mapToDouble(e -> e.err_obs * e.err_obs)
                                            .sum() / count
                            );
                            // 计算err_repair的RMSE
                            double rmseRepair = Math.sqrt(
                                    group.stream()
                                            .mapToDouble(e -> e.err_repair * e.err_repair)
                                            .sum() / count
                            );
                            // 返回比值
                            return rmseRepair / rmseObs;
                        }
                ));
    }
}

// 用于记录误差信息
class ErrorInfo {
    public long interval;       // 对于单向修复，是与前一个最近标记点时间戳之差，对于双向修复，是与最近标记点的时间戳之差。
    public double err_obs;
    public double err_repair;

    ErrorInfo(long interval, double err_obs, double err_repair) {
        this.interval = interval;
        this.err_obs = err_obs;
        this.err_repair = err_repair;
    }
}
