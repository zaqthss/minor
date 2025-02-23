package com.MINOR.experiment;

import com.MINOR.Utils.Constants;
import com.MINOR.Utils.MyFileUtils;
import com.MINOR.core.*;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimeSeries;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class Application {
    private enum Model {
        MINOR_B("MINOR-B"),
        MINOR_U("MINOR-U"),
        VARX("VARX"),
        IMR("IMR"),
        MTCSC("MTCSC"),
        Akane("Akane");

        private final String name;

        Model(String name) {
            this.name = name;
        }

        public String toString() {
            return this.name;
        }
    }

    public static void main(String[] args) {
        try {
            run("Car", "Classification", "TRAIN");
            run("Car", "Classification", "TEST");
            run("Lightning2", "Classification", "TRAIN");
            run("Lightning2", "Classification", "TEST");
            run("DistalPhalanxTW", "Cluster", "TRAIN");
            run("EOGVerticalSignal", "Cluster", "TRAIN");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void run(String dataset, String taskType, String suffix) throws IOException {
        Long[] seeds = Constants.seedsApp;
        int seedCnt = seeds.length;
        String cleanDataFile = Constants.rawPathPrefix + taskType + '/' + dataset + '/' + dataset + '_' + suffix + ".csv";
        List<List<Double>> cleanData = MyFileUtils.readAppData(cleanDataFile);
        String DirtyDataFile, LabelFile, outputFile;
        for (int i = 0; i < seedCnt; i++) {
            DirtyDataFile = Constants.appPathPrefix + "dirty/" + taskType + '/' + dataset + '/' + i + '_' + dataset + "_Dirty_" + suffix + ".csv";
            LabelFile = Constants.appPathPrefix + "dirty/" + taskType + '/' + dataset + '/' + "label_" + i + '_' + dataset + '_' + suffix + ".csv";
            Pair<List<Long>, List<List<Double>>> dirtyDara = readCsvFile(DirtyDataFile);
            List<List<Integer>> labels = MyFileUtils.readAppLabel(LabelFile);
            for (Model model : Model.values()) {
                List<List<Double>> repairedData = repair(dirtyDara.first, dirtyDara.second, cleanData, labels, model);
                outputFile = Constants.appPathPrefix + "repaired/" + taskType + '/' + dataset + '/' + i + '_' + dataset + '_' + model + '_' + suffix + ".csv";
                MyFileUtils.exportAppData(outputFile, dirtyDara.first, repairedData);
            }
        }
        System.out.println("[INFO]repair on " + dataset + '_' + suffix + " completed.");
    }

    private static List<List<Double>> repair(List<Long> timestamps, List<List<Double>> datas, List<List<Double>> cleanDatas, List<List<Integer>> labelList, Model model) {
        List<List<Double>> repairs = new ArrayList<>();
        int colNum = datas.size();
        for (int col = 0; col < colNum; col++) {
            List<Double> obs = datas.get(col);
            List<Double> truth = cleanDatas.get(col);
            List<Integer> lbs = labelList.get(col);
            TimeSeries ts = TimeSeries.getFromList(timestamps, obs, truth, lbs);
            List<Double> repair = new ArrayList<>();
            double max = Collections.max(obs);
            double min = Collections.min(obs);
            double mu = obs.stream().mapToDouble(Double::doubleValue).average().orElse(0);
            double thres = mu;
            int n = obs.size();
            int iterationLimit = 1_000;
            switch (model) {
                case MINOR_B: {
                    MINOR_B minorB = new MINOR_B(ts, 1, thres, iterationLimit);
                    TimeSeries repairedTs = minorB.runApp();
                    repair = repairedTs.toList();
                    break;
                }
                case MINOR_U: {
                    MINOR_U minorU = new MINOR_U(ts, 1, thres, iterationLimit);
                    TimeSeries repairedTs = minorU.runApp();
                    repair = repairedTs.toList();
                    break;
                }
                case VARX: {
                    VARX varx = new VARX(ts, 1, thres);
                    TimeSeries repairedTs = varx.runApp();
                    repair = repairedTs.toList();
                    break;
                }
                case IMR: {
                    UniTimeSeries uts = UniTimeSeries.getFromList(timestamps, obs, truth, lbs);
                    IMR imr = new IMR(uts, 1, thres, iterationLimit);
                    UniTimeSeries repairedTs = imr.runApp();
                    repair = repairedTs.toList();
                    break;
                }
                case MTCSC: {
                    double[] s = calculateMaxMinSpeed(obs, timestamps);
                    UniTimeSeries uts = UniTimeSeries.getFromList(timestamps, obs, truth, lbs);
                    MTCSC_Uni mtcsc = new MTCSC_Uni(uts, s[0], s[1], Math.min(n / 20, 100));
                    UniTimeSeries repairedTs = mtcsc.mainScreen();
                    repair = repairedTs.toList();
                    break;
                }
                case Akane: {
                    UniTimeSeries uts = UniTimeSeries.getFromList(timestamps, obs, truth, lbs);
                    AkaneHeuristic akaneheuristicauto = new AkaneHeuristic(uts, 3);
                    akaneheuristicauto.mainAkaneHeuristicAuto();
                    UniTimeSeries repairedTs = akaneheuristicauto.getResultTimeSeries();
                    repair = repairedTs.toList();
                    break;
                }
                default: {
                    System.out.println("[ERROR]no models are specified");
                }
            }
            repairs.add(repair);
        }
        return repairs;
    }

    private static double[] calculateMaxMinSpeed(List<Double> data, List<Long> timestamps) {
        if (data.size() != timestamps.size() || data.size() <= 1) {
            throw new IllegalArgumentException
                    ("[WARNING]The data list and the timestamp list must be of the same length and contain at least two elements.");
        }

        double maxSpeed = Double.NEGATIVE_INFINITY;
        double minSpeed = Double.POSITIVE_INFINITY;

        for (int i = 0; i < data.size() - 1; i++) {
            double deltaData = data.get(i + 1) - data.get(i);
            long deltaTime = timestamps.get(i + 1) - timestamps.get(i);

            if (deltaTime == 0) {
                throw new IllegalArgumentException
                        ("[WARNING]The timestamps are the same, resulting in a time difference of zero, which is not allowed.");
            }

            double speed = deltaData / deltaTime;

            if (speed > maxSpeed) {
                maxSpeed = speed;
            }
            if (speed < minSpeed) {
                minSpeed = speed;
            }
        }

        return new double[]{maxSpeed, minSpeed};
    }

    private static Pair<List<Long>, List<List<Double>>> readCsvFile(String csvFile) throws IOException {
        List<Long> timestamps = new ArrayList<>();
        List<List<Double>> datas = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                // First column as timestamps
                timestamps.add(Long.parseLong(values[0]));
                // Remaining columns as data
                for (int i = 1; i < values.length; i++) {
                    if (datas.size() < i) {
                        datas.add(new ArrayList<>());
                    }
                    datas.get(i - 1).add(Double.parseDouble(values[i]));
                }
            }
        }

        return new Pair<>(timestamps, datas);
    }

    private static class Pair<K, V> {
        public final K first;
        public final V second;

        public Pair(K first, V second) {
            this.first = first;
            this.second = second;
        }
    }
}
