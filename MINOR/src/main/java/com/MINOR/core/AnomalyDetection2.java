package com.MINOR.core;

import com.MINOR.entity.TimeSeries;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

@Deprecated
public class AnomalyDetection2 {

    public static List<Integer> detectAbnormalDimensions(TimeSeries data, double correlationThreshold, double dropThreshold) {
        int dimensions = data.size();
        double[][] correlationMatrix = new double[dimensions][dimensions];
        List<Integer> abnormalDimensions = new ArrayList<>();

        for (int i = 0; i < dimensions; i++) {
            int correlated = 0;
            int declined = 0;

            for (int j = 0; j < dimensions; j++) {
                if (i != j) {
                    // Calculate Pearson correlation between dimension i and j
                    double correlation = calculateCorrelation(data.getP(i).getValObs().data, data.getP(j).getValObs().data);
                    correlationMatrix[i][j] = correlation;

                    if (correlation >= correlationThreshold) {
                        correlated++;

                        // Check for significant drop in correlation
                        if (correlation < dropThreshold) {
                            declined++;
                        }
                    }
                }
            }

            // Mark dimension as abnormal if more than half of its correlations declined
            if (correlated > 0 && ((double) declined / correlated) > 0.5) {
                abnormalDimensions.add(i);
            }
        }
        return abnormalDimensions;
    }

    private static double calculateCorrelation(double[] dim1, double[] dim2) {
        double mean1 = Arrays.stream(dim1).average().orElse(0);
        double mean2 = Arrays.stream(dim2).average().orElse(0);

        double numerator = 0.0, denom1 = 0.0, denom2 = 0.0;
        for (int i = 0; i < dim1.length; i++) {
            double diff1 = dim1[i] - mean1;
            double diff2 = dim2[i] - mean2;
            numerator += diff1 * diff2;
            denom1 += diff1 * diff1;
            denom2 += diff2 * diff2;
        }
        return numerator / Math.sqrt(denom1 * denom2);
    }
}
