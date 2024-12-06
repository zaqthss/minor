package com.MINOR.pre;

import com.MINOR.Utils.Constants;

public class HandleGPS {
    public static void main(String[] args) {
        double[] lrs = Constants.lrs;
        String name = "GPS";
        String inputFPath = Constants.rawPathPrefix + "GPS.csv";
        int caseCnt = Constants.lrs.length;
        for (int i = 0; i < caseCnt; i++) {
            for (int seed = 0; seed < Constants.seeds.length; seed++) {
                String outputFPath = Constants.dirtyPathPrefix + "GPS/GPS_lr_" + lrs[i] + '_' + seed + ".csv";
                GPSHandler handler = new GPSHandler(inputFPath, outputFPath, lrs[i], Constants.seeds[seed]);
                handler.run();
            }
        }
        System.out.println("[DONE]label preprocess for " + name + " done");
    }
}
