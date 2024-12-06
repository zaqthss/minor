package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.TimeSeries;
import org.ejml.dense.row.MatrixFeatures_DDRM;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;


/**
 * 异常检测接口，实现该接口即可继承异常检测方法
 */
@Deprecated
public interface AnomalyDetector {
    /**
     * 对时序数据进行异常检测，将检测出的正常点进行标记
     * 该方法是calculateSpeedThreshold与performAnomalyDetection的封装
     *
     * @param ts 需要进行异常检测的时序数据
     */
    default void detectAnomalies(TimeSeries ts) {
        performAnomalyDetection(ts, calculateSpeedThreshold(ts));
    }

    /**
     * 计算速度约束
     *
     * @param ts 时序数据
     * @return 速度约束
     */
    default double calculateSpeedThreshold(TimeSeries ts) {
        List<Double> speeds = new ArrayList<>();
        int dim = ts.getDimension();
        for (int i = 1; i < ts.size(); i++) {
            double distance = MyMathUtils.calL2Norm(ts.getObsAt(i), ts.getObsAt(i - 1), dim);
            double timeDiff = Math.abs(ts.getTimestampAt(i) - ts.getTimestampAt(i - 1));
            double speed = distance / timeDiff;
            speeds.add(speed);
        }

        double meanSpeed = speeds.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double variance = speeds.stream().mapToDouble(speed -> pow(speed - meanSpeed, 2)).average().orElse(0.0);
        double stdDev = sqrt(variance);

        return (meanSpeed + 2 * stdDev);
    }

    /**
     * 对时序数据进行异常检测，将检测出的正常点进行标记
     *
     * @param ts   需要进行检测的时序数据
     * @param Smax 速度约束
     */
    default void performAnomalyDetection(TimeSeries ts, double Smax) {
        int dim = ts.getDimension();
        int n = ts.size();
        int[] anomaly = new int[n];
        int[] normal = new int[n];
        for (int i = 0; i < n; i++) {
            anomaly[i] = i - 1;
            normal[i] = -1;
        }
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < j; i++) {
                double speed = MyMathUtils.calL2Norm(ts.getObsAt(i), ts.getObsAt(j), dim) / (ts.getTimestampAt(j) - ts.getTimestampAt(i));
                if (speed <= Smax && anomaly[j] > anomaly[i] + (j - i - 1)) {
                    anomaly[j] = anomaly[i] + (j - i - 1);
                    normal[j] = i;
                }
            }
        }

        int index = anomaly[0] + (n - 1);
        int res = 0;

        for (int j = 1; j < n; j++) {
            if (anomaly[j] + (n - j) < index) {
                index = anomaly[j] + (n - j);
                res = j;
            }
        }
        for (int x : normal) {
            if (x >= 0) {
                ts.setDetectedAsNormal(x, true);
            }
        }
        ts.applyAD();
    }

    /**
     * 衡量异常检测的性能
     *
     * @param ts    进行异常检测后的时序数据
     */
    default void ADEvaluate(TimeSeries ts) {
        int dim = ts.getDimension();
        int n = ts.size();
        int TP = 0;
        int FP = 0;
        int TN = 0;
        int FN = 0;
        double accuracy, precision, recall, F1;
        for (int i = 0; i < n; i++) {
            if (ts.isDetectedAsNormal(i)) {
                if (MatrixFeatures_DDRM.isEquals(ts.getObsAt(i), ts.getTruthAt(i))) {
                    TN++;
                } else {
                    FN++;
                }
            } else {
                if (MatrixFeatures_DDRM.isEquals(ts.getObsAt(i), ts.getTruthAt(i))) {
                    FP++;
                } else {
                    TP++;
                }
            }
        }
        accuracy = (double) (TP + TN) / n;
        precision = (TP + FP) != 0 ? (double) TP / (TP + FP) : 0;
        recall = (TP + FN) != 0 ? (double) TP / (TP + FN) : 0;
        F1 = (precision + recall) != 0 ? 2 * (precision * recall) / (precision + recall) : 0;
        System.out.println("[INFO]:\n\tAccuracy=" + accuracy + "\n\tprecision=" + precision + "\n\trecall=" + recall + "\n\tF1=" + F1);
    }
}
