package com.MINOR.experiment;

import com.MINOR.Utils.Constants;
import com.MINOR.Utils.MyFileUtils;
import com.MINOR.core.*;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimeSeries;

import java.util.List;

/**
 * repaired time series data of Fig.4
 */
public class CaseStudy_Figure {
    public static void main(String[] args) {
        final double lr0 = Constants.lr0;
        final int seed = 0;

        final String dataFPath = Constants.dirtyPathPrefix + "GPS/GPS_lr_" + lr0 + '_' + seed + ".csv";

        TimeSeries ts = new TimeSeries();
        int dimension = 2;

        int p = 1;
        double threshold = 0.05;
        int iterationLimit = 20000;
        int T = 100;
        double SMax = 1.4325;

        // MINOR-B
        String resultFilename = Constants.resultPathPrefix + "CaseStudy/MINOR-B.csv";
        MyFileUtils.csv2MTS(ts, dimension, dataFPath);
        MINOR_B minorB = new MINOR_B(ts, p, threshold, iterationLimit);
        minorB.run();
        minorB.exportTimeSeries(resultFilename);
        // MINOR-U
        resultFilename = Constants.resultPathPrefix + "CaseStudy/MINOR-U.csv";
        MyFileUtils.csv2MTS(ts, dimension, dataFPath);
        MINOR_U minorU = new MINOR_U(ts, p, threshold, iterationLimit);
        minorU.run();
        minorU.exportTimeSeries(resultFilename);
        // IMR
        resultFilename = Constants.resultPathPrefix + "CaseStudy/IMR.csv";
        MyFileUtils.csv2MTS(ts, dimension, dataFPath);
        List<UniTimeSeries> unis = ts.toUniTimeSeries();
        for (int j = 0; j < dimension; j++) {
            IMR imr = new IMR(unis.get(j), 1, threshold, iterationLimit);
            imr.run();
        }
        TimeSeries aggTs = TimeSeries.aggregate(unis);
        MyFileUtils.MTSExporter(resultFilename, aggTs);
        // MTCSC-A
        resultFilename = Constants.resultPathPrefix + "CaseStudy/MTCSC-A.csv";
        MyFileUtils.csv2MTS(ts, dimension, dataFPath);
        MTCSC_A mtcsc_a = new MTCSC_A(ts, SMax, T, 0.10, 0.75, 125, 0.75, 5);
        mtcsc_a.mainScreen();
        MyFileUtils.MTSExporter(resultFilename, ts);

    }
}
