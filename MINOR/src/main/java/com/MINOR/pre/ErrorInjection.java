package com.MINOR.pre;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.enums.ErrorInjectorType;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class ErrorInjection {
    private double mu, sigma;
    private int dataNum;
    private int dim;        // dimension of the data

    private double lRate; // 标记数据占比
    private int eSections, eLength;  // 生成错误的起点数量,连续错误的长度
    private List<Integer> startPts;


    private List<List<Double>> data;
    private List<Integer> timestamps;
    private double[] mean;      // 原始数据均值
    private double[][] covMatrix;   // 原始数据协方差
    private double[][] corMatrix;   //  原始数据相关系数矩阵
    private ErrorInjectorType et;

    private long seed;
    private String outputFPath;

    /**
     * @param _data        raw data
     * @param _timestamps  timestamp list
     * @param _mean        mean of the raw data
     * @param _covMatrix   covariance matrix of the raw data
     * @param _corMatrix   correlation matrix of the raw data
     * @param _mu          mu of the Gaussian Noise
     * @param _sigma       sigma of the Gaussian Noise
     * @param _lRate       label rate
     * @param _eLength     error length
     * @param _eSections   the number of error sections
     * @param _seed        random seed
     * @param _et          type of the error that intend to be injected, SHIFT or INNOVATIONAL
     * @param _outputFPath the file path of the data that injected errors
     */
    public ErrorInjection(List<List<Double>> _data, List<Integer> _timestamps, double[] _mean, double[][] _covMatrix, double[][] _corMatrix,
                          double _mu, double _sigma, double _lRate, int _eLength, int _eSections, ErrorInjectorType _et,
                          long _seed, String _outputFPath) {
        data = _data;
        timestamps = _timestamps;
        mean = _mean;
        covMatrix = _covMatrix;
        corMatrix = _corMatrix;
        mu = _mu;
        sigma = _sigma;
        lRate = _lRate;
        eLength = _eLength;
        eSections = _eSections;
        seed = _seed;
        et = _et;
        outputFPath = _outputFPath;
    }

    public void run() {
        double[] meanScaled = MyMathUtils.scaleMeanByMu(this.mean, this.mu);
        double[][] cov = MyMathUtils.generateCovByCor(this.corMatrix, Math.pow(this.sigma, 2));
        dataNum = data.size();
        dim = data.get(0).size();
        startPts = new ArrayList<>();

        Random random = new Random();
        random.setSeed(seed);

        List<List<Double>> dirtyData = new ArrayList<>();  // 加错后的数据
        List<Double> rowCopy;
        Set<Integer> lbSet = new HashSet<>();
        // 深拷贝数据
        for (int i = 0; i < dataNum; i++) {
            rowCopy = new ArrayList<>(data.get(i));
            dirtyData.add(rowCopy);
        }
        // 计算每个段的长度
        int sectionLength = dataNum / eSections;
        // 添加噪音
        for (int i = 0; i < eSections; i++) {
            // 随机选择一个起始点，但保证起始点距离段末尾足够远
            int start = random.nextInt(sectionLength - eLength - 2) + i * sectionLength;
            startPts.add(start);
            int markedPointsCount = 0;
            Set<Integer> markedIndices = new HashSet<>();
            if (eLength != 1) {
//                if (eLength <= 5) {
//                    // 只标记起点
//                    markedPointsCount = Math.max(1, (int) Math.ceil(lRate * eLength));
//                    markedIndices.add(start);
//                    lbSet.add(start);
//                } else {
//                    // 计算需要标记的点的数量，确保包括前两个点
//                    markedPointsCount = Math.max(2, (int) Math.ceil(lRate * eLength));
//                    // 默认标记起始点和起始点+1
//                    markedIndices.add(start);
//                    markedIndices.add(start + 1);
//                    lbSet.add(start);
//                    lbSet.add(start + 1);
//                }
                // 一定标记起点
                markedPointsCount = Math.max(1, (int) Math.floor(lRate * eLength));
                markedIndices.add(start);
                lbSet.add(start);
                // 选择随机点进行标记，保证不重复，且数量满足需求
                while (markedIndices.size() < markedPointsCount) {
                    int markIndex = start + random.nextInt(eLength);
                    markedIndices.add(markIndex);
                    lbSet.add(markIndex);
                }
            }
            int sign = random.nextInt(2) * 2 - 1;// 噪音的符号正负
            int type = random.nextInt(3);   // subtypes of INNOVATIONAL error
            // 对每个段添加噪音
            for (int j = 0; j < eLength; j++) {
                int index = start + j;
                List<Double> errs = MyMathUtils.generateGaussianNoise(meanScaled, cov, dim, random);
                if (et == ErrorInjectorType.INNOVATIONAL) {
                    double gradient = 1;
                    switch (type) {
                        case 0: {
                            // UP
                            gradient = (double) j / eLength;
                            break;
                        }
                        case 1: {
                            // DOWN
                            gradient = 1 - (double) j / eLength;
                            break;
                        }
                        case 2: {
                            // MID
                            int mid = eLength / 2;
                            gradient = 1 - (double) Math.abs(j - mid) / mid;
                            break;
                        }
                    }
                    for (int k = 0; k < errs.size(); k++) {
                        errs.set(k, errs.get(k) * gradient);
                    }
                }
                for (int m = 0; m < dim; m++) {
                    dirtyData.get(index).set(m, dirtyData.get(index).get(m) + sign * errs.get(m));
                }
            }
        }
        if (eLength == 1){
            // 对于spike errors, startPts就是所有脏值点
            int markedPointsCount = Math.max(1, (int) Math.floor(lRate * eSections));
            Collections.shuffle(startPts);
            for (int i = 0; i < markedPointsCount; i++) {
                lbSet.add(startPts.get(i));
            }
        }


        CSVFormat csvFormat = CSVFormat.DEFAULT.builder()
                .setSkipHeaderRecord(false)
                .build();
        try (
                BufferedWriter dataWriter = Files.newBufferedWriter(Paths.get(outputFPath));
                CSVPrinter csvPrinter = new CSVPrinter(dataWriter, csvFormat)
        ) {
            // 写数据，按照时间戳，加错后数据，原始数据，标签的顺序
            for (int i = 0; i < dataNum; i++) {
                List<Object> records = getObjects(i, dirtyData, lbSet);
                csvPrinter.printRecord(records);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private List<Object> getObjects(int i, List<List<Double>> dirtyData, Set<Integer> lbSet) {
        List<Object> records = new ArrayList<>();
        records.add(timestamps.get(i));
        for (int j = 0; j < dim; j++) {
            records.add(dirtyData.get(i).get(j));
        }
        for (int j = 0; j < dim; j++) {
            records.add(data.get(i).get(j));
        }
        if (lbSet.contains(i)) {
            records.add(1);
        } else {
            records.add(0);
        }
        return records;
    }

    public void setMu(double mu) {
        this.mu = mu;
    }

    public void setSigma(double sigma) {
        this.sigma = sigma;
    }

    public void seteLength(int eLength) {
        this.eLength = eLength;
    }

    public void setlRate(double lRate) {
        this.lRate = lRate;
    }
}
