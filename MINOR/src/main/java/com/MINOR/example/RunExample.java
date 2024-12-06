package com.MINOR.example;

import com.MINOR.Utils.Constants;
import com.MINOR.Utils.MyFileUtils;
import com.MINOR.core.*;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimeSeries;

import java.util.List;

public class RunExample {
    public static void main(String[] args) {
        String csvFilePath = Constants.dirtyPathPrefix + "example/example.csv";
        TimeSeries ts = new TimeSeries();
        int dimension = 2;
        MyFileUtils.csv2MTS(ts, dimension, csvFilePath);

        int p = 1;
        double threshold = 0.02;
        int iterationLimit = 10000;

        // MINOR-B
        String resultFilename = Constants.resultPathPrefix + "example/example_MINOR_B.csv";
        MINOR_B minorbi = new MINOR_B(ts, p, threshold, iterationLimit);
        minorbi.run();
        minorbi.exportTimeSeries(resultFilename);
        // MINOR-U
        MyFileUtils.csv2MTS(ts, dimension, csvFilePath);
        resultFilename = Constants.resultPathPrefix + "example/example_MINOR_U.csv";
        MINOR_U minoru = new MINOR_U(ts, p, threshold, iterationLimit);
        minoru.runWithoutCandiIC();
        minoru.exportTimeSeries(resultFilename);
        // VAR
        MyFileUtils.csv2MTS(ts, dimension, csvFilePath);
        resultFilename = Constants.resultPathPrefix + "example/example_VAR.csv";
        VAR var = new VAR(ts, p, threshold);
        var.run();
        var.exportTimeSeries(resultFilename);
        // VARX
        MyFileUtils.csv2MTS(ts, dimension, csvFilePath);
        resultFilename = Constants.resultPathPrefix + "example/example_VARX.csv";
        VARX varx = new VARX(ts, p, threshold);
        varx.run();
        varx.exportTimeSeries(resultFilename);
        // IMR
        resultFilename = Constants.resultPathPrefix + "example/example_IMR.csv";
        MyFileUtils.csv2MTS(ts, dimension, csvFilePath);
        List<UniTimeSeries> unis = ts.toUniTimeSeries();
        for (int j = 0; j < 2; j++) {
            IMR imr = new IMR(unis.get(j), p, threshold, iterationLimit);
            imr.run();
        }
        TimeSeries aggTs = TimeSeries.aggregate(unis);
        MyFileUtils.MTSExporter(resultFilename, aggTs);
    }
}
