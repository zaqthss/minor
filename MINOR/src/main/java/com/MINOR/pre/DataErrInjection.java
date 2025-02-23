package com.MINOR.pre;

import com.MINOR.Utils.Constants;
import com.MINOR.Utils.MyFileUtils;
import com.MINOR.Utils.MyMathUtils;
import com.MINOR.enums.ErrorInjectorType;

import java.util.ArrayList;
import java.util.List;

/**
 * 对数据集注入错误值
 */
public class DataErrInjection {
    private static double lr0 = 0.2;

    public static void main(String[] args) {
        injector("ILD", 9603, 2, 3, 0.1);
        injector("ECG", 9394, 32, 1, 0.05);
    }

    private static void injector(String name, int dataSize, int dim, double mu, double sigma) {
        // raw data
        String filePath = Constants.rawPathPrefix + name + ".csv";
        List<List<Double>> data = MyFileUtils.csv2List(filePath, 1, dim);
        List<Integer> timestamps = MyFileUtils.csv2IntList(filePath, 0);
        inject(data, timestamps, dataSize, mu, sigma, name);
        List<List<Double>> dataFull;
        int[] sizes;
        if (name.equals("ILD")) {
            filePath = Constants.rawPathPrefix + "ILD48k.csv";
            sizes = Constants.sizesILD;
        } else {
            filePath = Constants.rawPathPrefix + "ECG93k.csv";
            sizes = Constants.sizesECG;
        }
        timestamps = MyFileUtils.csv2IntList(filePath, 0);
        dataFull = MyFileUtils.csv2List(filePath, 1, dim);
        truncateAndInject(dataFull, timestamps, sizes, mu, sigma, name);
        System.out.println("[DONE]Error injection for " + name + " done");
    }

    /**
     * error rate and error length experiment
     *
     * @param data
     * @param timestamps
     * @param dataSize
     * @param mu
     * @param sigma
     * @param name
     */
    private static void inject(List<List<Double>> data, List<Integer> timestamps, int dataSize, double mu, double sigma, String name) {
        double[] mean = MyMathUtils.calculateMean(data);
        double[][] cov = MyMathUtils.calculateCovariance(data, mean);
        double[][] cor = MyMathUtils.calculateCorrelation(cov);
        double er0 = Constants.er0;
        int el0 = Constants.el0;
        double[] ers = Constants.ers;
        int[] els = Constants.els;
        for (int seed = 0; seed < Constants.seeds.length; seed++) {
            // 改变错误率
            for (int i = 0; i < 10; i++) {
                // 有多少段错误
                int sections = (int) Math.round((double) dataSize * ers[i] / el0);
                String outFP = Constants.dirtyPathPrefix + "\\" + name + "\\" + name + "_SHIFT_er_" + ers[i] + '_' + seed + ".csv";
                ErrorInjection err = new ErrorInjection(data, timestamps, mean, cov, cor, mu, sigma,
                        lr0, el0, sections, ErrorInjectorType.SHIFT,
                        Constants.seeds[seed], outFP);
                err.run();
                if(name.equals("ILD")) {
                    // error pattern experiment will only be conducted on the ILD dataset
                    String outFP2 = Constants.dirtyPathPrefix + "\\" + name + "\\" + name + "_INNO_er_" + ers[i] + '_' + seed + ".csv";
                    ErrorInjection err2 = new ErrorInjection(data, timestamps, mean, cov, cor, mu, sigma,
                            lr0, el0, sections, ErrorInjectorType.INNOVATIONAL,
                            Constants.seeds[seed], outFP2);
                    err2.run();
                }
            }
            // 改变长度
            for (int i = 0; i < 10; i++) {
                int sections = (int) Math.round((double) dataSize * er0 / els[i]);
                String outFP = Constants.dirtyPathPrefix + "\\" + name + "\\" + name + "_SHIFT_el_" + els[i] + '_' + seed + ".csv";
                ErrorInjection err = new ErrorInjection(data, timestamps, mean, cov, cor, mu, sigma,
                        lr0, els[i], sections, ErrorInjectorType.SHIFT,
                        Constants.seeds[seed], outFP);
                err.run();
            }
        }

    }

    /**
     * data size experiment
     *
     * @param data
     * @param timestamps
     * @param sizes
     * @param mu
     * @param sigma
     * @param name
     */
    private static void truncateAndInject(List<List<Double>> data, List<Integer> timestamps, int[] sizes, double mu, double sigma, String name) {
        for (int seed = 0; seed < Constants.seeds.length; seed++) {
            for (int i = 0; i < sizes.length; i++) {
                List<Double> rowCopy;
                List<List<Double>> subList = new ArrayList<>();
                for (int j = 0; j < sizes[i]; j++) {
                    rowCopy = new ArrayList<>(data.get(j));
                    subList.add(rowCopy);
                }
                double[] mean = MyMathUtils.calculateMean(subList);
                double[][] cov = MyMathUtils.calculateCovariance(subList, mean);
                double[][] cor = MyMathUtils.calculateCorrelation(cov);
                double er0 = Constants.er0;
                int el0 = Constants.el0;

                int sections = (int) Math.round((double) sizes[i] * er0 / el0);
                String outFP = Constants.dirtyPathPrefix + "\\" + name + "\\" + name + "_SHIFT_size_" + sizes[i] + '_' + seed + ".csv";
                ErrorInjection err = new ErrorInjection(subList, timestamps, mean, cov, cor, mu, sigma,
                        lr0, el0, sections, ErrorInjectorType.SHIFT,
                        Constants.seeds[seed], outFP);
                err.run();
            }
        }
    }
}
