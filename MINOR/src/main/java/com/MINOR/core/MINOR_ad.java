package com.MINOR.core;

import com.MINOR.Utils.MyFileUtils;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

@Deprecated
public class MINOR_ad extends MINOR_base {
    private double Smax;
    protected List<Double> speeds;

    public MINOR_ad(TimeSeries ts, int p, double threshold, long iterationLimit) {
        this.modelName = "MONOR-AD";

        this.ts = ts;
        this.p = p;
        this.threshold = threshold;
        this.iterationLimit = iterationLimit;

        r = -1;
        this.dim = ts.getDimension();
        this.tsLength = ts.size();
        Z = new DMatrixRMaj(tsLength - p, dim * p);
        V = new DMatrixRMaj(tsLength - p, dim);
        this.phi = new DMatrixRMaj(p * dim, dim);
    }

    public void calcsmax() {
        for (int i = 1; i < ts.size(); i++) {
            double distance = sqrt(pow(ts.s.get(i).getValObs().get(0) - ts.s.get(i - 1).getValObs().get(0), 2) +
                    pow(ts.s.get(i).getValObs().get(1) - ts.s.get(i - 1).getValObs().get(1), 2));
            double timeDiff = ts.s.get(i).getTimestamp() - ts.s.get(i - 1).getTimestamp();
            if (timeDiff > 0) {
                double speed = distance / timeDiff;
                speeds.add(speed);
            }
        }

        double meanSpeed = speeds.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double variance = speeds.stream().mapToDouble(speed -> pow(speed - meanSpeed, 2)).average().orElse(0.0);
        double stdDev = sqrt(variance);

        this.Smax = meanSpeed + 2 * stdDev;
    }

    public void anomalyDetection() {
        int n = ts.size();
        int[] anomaly = new int[n];
        int[] normal = new int[n];
        for (int i = 0; i < n; i++) {
            anomaly[i] = i - 1;
            normal[i] = -1;
        }
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < j; i++) {
                double speed = sqrt(pow((ts.getObsAt(j).get(0) - ts.getObsAt(i).get(0)), 2) + pow((ts.getObsAt(j).get(1) - ts.getObsAt(i).get(1)), 2)) / (ts.getTimestampAt(j) - ts.getTimestampAt(i));
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
        for (int i = 0; i < n; i++) {
            if (i == normal[res]) {
                ts.s.get(i).setLabel(true);
            }
        }

    }

    public TimeSeries run(String resultFilename, String phiFilename) {
        this.speeds = new ArrayList<>();
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        preprocess();
        MatrixPruning();
        int k = 1;
        this.isConverge = false;
        calcsmax();
        anomalyDetection();
        while (true) {
            estimateIC();
            candidate();
            evaluate();
            if (isConverge || k >= iterationLimit) {
                System.out.println("[INFO]Done " + this.modelName + "_ad(" + p + ") stopped at " + k + "th iteration.");
                break;
            }
            k++;
        }
        MyFileUtils.MTSExporter(resultFilename, ts);
        MyFileUtils.PhiExporter(phiFilename, this.phi);
        return ts;
    }

    /**
     * 只进行异常检测，将正常点的label设为true
     *
     * @return 标记过的timeSeries
     */
    public TimeSeries ad() {
        this.speeds = new ArrayList<>();
        A = new DMatrixRMaj();
        B = new DMatrixRMaj();
        preprocess();
        MatrixPruning();
        this.isConverge = false;
        calcsmax();
        anomalyDetection();
        return ts;
    }
}