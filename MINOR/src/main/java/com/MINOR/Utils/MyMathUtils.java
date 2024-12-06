package com.MINOR.Utils;

import com.MINOR.entity.TimePoint;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimeSeries;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.math3.linear.*;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.MatrixFeatures_DDRM;

import java.io.*;
import java.util.*;

public class MyMathUtils {
    /**
     * 计算两个点之间的L2范数
     *
     * @param point1    列向量
     * @param point2    列向量
     * @param dimension 向量的维度
     * @return L2范数
     */
    public static double calL2Norm(DMatrixRMaj point1, DMatrixRMaj point2, int dimension) {
        if (dimension <= 0 || dimension > point1.numRows || dimension > point2.numRows) {
            throw new IllegalArgumentException("维度必须在有效范围内");
        }
        double sum = 0.0;
        for (int i = 0; i < dimension; i++) {
            double diff = point1.get(i, 0) - point2.get(i, 0);
            sum += Math.pow(diff, 2);
        }
        return Math.sqrt(sum);
    }

    /**
     * 计算多维时序数据两个数据点的修复值之间的速度，不需要考虑两个参数的先后
     */
    public static double calSpeed(TimePoint point1, TimePoint point2) {
        return calL2Norm(point1.getValRepaired(), point2.getValRepaired(), point1.getDimension()) / Math.abs(point1.getTimestamp() - point2.getTimestamp());
    }

    /**
     * 计算两个向量夹角的余弦值
     *
     * @param vec1      当前修复值-上一次修复值
     * @param vec2      candidate-当前修复值
     * @param dimension 时序数据的维度
     * @return vec1与vec2的夹角余弦值
     */
    public static double cosineSimilarity(DMatrixRMaj vec1, DMatrixRMaj vec2, int dimension) {
        if (dimension <= 0 || dimension > vec1.numRows || dimension > vec2.numRows || vec1.numCols > 1 || vec2.numCols > 1) {
            throw new IllegalArgumentException("维度必须在有效范围内");
        }
        double l2n1 = 0.0;
        double l2n2 = 0.0;
        double dotProduct;
        for (int i = 0; i < dimension; i++) {
            l2n1 += Math.pow(vec1.get(i, 0), 2);
            l2n2 += Math.pow(vec2.get(i, 0), 2);
        }
        l2n1 = Math.sqrt(l2n1);
        l2n2 = Math.sqrt(l2n2);
        if (l2n1 * l2n2 == 0) {
            return 0;
        }
        // 该函数仅允许列向量
        dotProduct = CommonOps_DDRM.dot(vec1, vec2);
        return dotProduct / (l2n1 * l2n2);
    }

    /**
     * vec1=repair-lastModify, vec2=candidate-repair,计算vec1和vec2的夹角余弦值
     *
     * @param lastModify 上一次的值，当且仅当在被修复过一次后这个值才会变更
     * @param repair     当前修复值
     * @param candidate  candidate
     * @param dimension  时序数据维度
     * @return vec1与vec2的夹角余弦值
     */
    public static double cosineSimilarityOfRepair(DMatrixRMaj lastModify, DMatrixRMaj repair, DMatrixRMaj candidate, int dimension) {
        if (dimension <= 0 || dimension > repair.numRows || repair.numCols > 1) {
            throw new IllegalArgumentException("维度必须在有效范围内");
        }
        DMatrixRMaj vec1 = new DMatrixRMaj(repair.numRows, repair.numCols);
        CommonOps_DDRM.subtract(repair, lastModify, vec1);
        DMatrixRMaj vec2 = new DMatrixRMaj(repair.numRows, repair.numCols);
        CommonOps_DDRM.subtract(candidate, repair, vec2);

        return cosineSimilarity(vec1, vec2, dimension);
    }

    /**
     * 获得时序数据某一个点本次修改和上一次修改向量的夹角
     *
     * @param ts    时序数据
     * @param index 序号
     * @return 两次修改向量的夹角余弦值
     */
    public static double cosineSimilarityOfRepair(TimeSeries ts, int index) {
        return cosineSimilarityOfRepair(ts.getLastModifyAt(index), ts.getRepairAt(index), ts.getCandidateAt(index), ts.getDimension());
    }

    /**
     * for uni-variate version
     */
    public static boolean isRoundTrip(UniTimeSeries ts, int index) {
        return (ts.getRepairAt(index) - ts.getLastModifyAt(index)) * (ts.getCandidateAt(index) - ts.getRepairAt(index)) < 0;
    }

    /**
     * 将两个一维序列合成二维序列
     */
    @Deprecated
    public static void mergeListsToMatrices(List<Double> x1, List<Double> x2, List<DMatrixRMaj> x) {
        if (x1.size() != x2.size()) {
            throw new IllegalArgumentException("Lists x1 and x2 must have the same size.");
        }
        x.clear();
        for (int i = 0; i < x1.size(); i++) {
            DMatrixRMaj matrix = new DMatrixRMaj(2, 1);
            matrix.set(0, 0, x1.get(i));
            matrix.set(1, 0, x2.get(i));
            x.add(matrix);
        }
    }

    /**
     * 删除矩阵Z全0的行，并删去V对应行
     */
    public static void MatrixPruning(DMatrixRMaj[] wrapper) {
        DMatrixRMaj Z = wrapper[0];
        DMatrixRMaj V = wrapper[1];
        int rows = Z.numRows;
        int cols = Z.numCols;
        int nonZeroRowCount = 0;
        boolean[] isZeroRow = new boolean[rows];
        // 遍历矩阵的每一行
        for (int i = 0; i < rows; i++) {
            boolean allZero = true;
            // 检查这一行的每一个元素
            for (int j = 0; j < cols; j++) {
                if (Z.get(i, j) != 0) {
                    allZero = false;
                    break;
                }
            }
            if (allZero) {
                isZeroRow[i] = true; // 标记这一行是全为0的行
            } else {
                nonZeroRowCount++;
            }
        }
        // 创建新的矩阵存储非零行
        DMatrixRMaj Z_pruned = new DMatrixRMaj(nonZeroRowCount, cols);
        DMatrixRMaj V_pruned = new DMatrixRMaj(nonZeroRowCount, V.numCols);
        int newRowIndex = 0;
        // 复制非零行到新矩阵
        for (int i = 0; i < rows; i++) {
            if (!isZeroRow[i]) {
                for (int j = 0; j < cols; j++) {
                    Z_pruned.set(newRowIndex, j, Z.get(i, j));
                }
                for (int k = 0; k < V.numCols; k++) {
                    V_pruned.set(newRowIndex, k, V.get(i, k));
                }
                newRowIndex++;
            }
        }
        wrapper[0] = Z_pruned;
        wrapper[1] = V_pruned;
    }

    /**
     * 以L2范数计算两个时序数据的欧式距离，以观察值为准
     */
    public static double eulerDistObs(TimeSeries ts1, TimeSeries ts2) {
        double dist = 0;
        int n = ts1.size();
        int m = ts1.getDimension();
        if (n != ts2.size()) {
            System.err.println("[ERROR]: the size of two time-series doesn't match.");
            System.exit(1);
        }
        for (int i = 0; i < n; i++) {
            dist += Math.pow(MyMathUtils.calL2Norm(ts1.getObsAt(i), ts2.getObsAt(i), m), 2);
        }
        dist = Math.sqrt(dist);
        return dist;
    }

    /**
     * 数字转布尔
     *
     * @param x csv文件中int形式的label
     * @return 1 true
     */
    public static boolean int2Bool(int x) {
        return x != 0;
    }

    /**
     * 布尔转数字，方便打印
     *
     * @param b label
     * @return true 1
     */
    public static int bool2Int(boolean b) {
        return b ? 1 : 0;
    }

    /**
     * 对异常检测算法进行评分，打印Accuracy, Precision, Recall, F1
     *
     * @param ts 进行异常检测后的时序数据
     */
    public static void ADEvaluate(TimeSeries ts) {
        int n = ts.size();
        int TP = 0;
        int FP = 0;
        int TN = 0;
        int FN = 0;
        double accuracy, precision, recall, F1;
        for (int i = 0; i < n; i++) {
            if (ts.getAnomalyAt(i)) {
                if (MatrixFeatures_DDRM.isEquals(ts.getObsAt(i), ts.getTruthAt(i))) {
                    FP++;
                } else {
                    TP++;
                }
            } else {
                if (MatrixFeatures_DDRM.isEquals(ts.getObsAt(i), ts.getTruthAt(i))) {
                    TN++;
                } else {
                    FN++;
                }
            }
        }
        accuracy = (double) (TP + TN) / n;
        precision = (TP + FP) != 0 ? (double) TP / (TP + FP) : 0;
        recall = (TP + FN) != 0 ? (double) TP / (TP + FN) : 0;
        F1 = (precision + recall) != 0 ? 2 * (precision * recall) / (precision + recall) : 0;
        System.out.println("[INFO]:\n\tAccuracy=" + accuracy + "\n\tprecision=" + precision + "\n\trecall=" + recall + "\n\tF1=" + F1);
    }

    /**
     * Generate n-dimensional Gaussian noise.
     *
     * @param mean       The mean vector
     * @param covariance The covariance matrix
     * @param dim        The dimension
     * @param rand       A random number generator with a specified seed
     * @return n-dimensional Gaussian noise
     */
    public static List<Double> generateGaussianNoise(double[] mean, double[][] covariance, int dim, Random rand) {
        RealMatrix covMatrix = new Array2DRowRealMatrix(covariance);
        // Cholesky分解
        CholeskyDecomposition decomposition = new CholeskyDecomposition(covMatrix);
        RealMatrix cholMatrix = decomposition.getL();
        // 生成标准正态分布的随机样本
        double[] standardNormal = new double[dim];
        for (int i = 0; i < dim; i++) {
            standardNormal[i] = rand.nextGaussian();
        }
        // 将标准正态分布的随机样本转换为具有指定协方差的分布
        RealVector sampleVector = new ArrayRealVector(standardNormal);
        RealVector transformedSample = cholMatrix.operate(sampleVector);
        // 将转换后的样本加上均值，得到最终的高斯噪声
        for (int i = 0; i < dim; i++) {
            transformedSample.setEntry(i, transformedSample.getEntry(i) + mean[i]);
        }
        // 将结果转换为List<Double>
        List<Double> gaussianNoise = new ArrayList<>();
        for (int i = 0; i < dim; i++) {
            gaussianNoise.add(transformedSample.getEntry(i));
        }
        return gaussianNoise;
    }

    /**
     * Calculate the mean of each dimension in a multidimensional dataset.
     */
    public static double[] calculateMean(List<List<Double>> data) {
        int n = data.size();
        int dim = data.get(0).size();
        double[] means = new double[dim];
        for (int i = 0; i < dim; i++) {
            means[0] = 0.0;
        }
        for (List<Double> row : data) {
            for (int i = 0; i < dim; i++) {
                if (i < row.size() && row.get(i) != null) {
                    means[i] = row.get(i) + means[i];
                }
            }
        }
        for (int i = 0; i < dim; i++) {
            means[i] /= n;
        }
        return means;
    }

    /**
     * Calculate the covariance matrix of a multidimensional dataset
     *
     * @param data the data
     * @param mean the sequence of means
     * @return the covariance matrix
     */
    public static double[][] calculateCovariance(List<List<Double>> data, double[] mean) {
        int n = data.size();
        int dim = data.get(0).size();
        double[][] covarianceMatrix = new double[dim][dim];
        for (List<Double> point : data) {
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                    covarianceMatrix[j][k] += (point.get(j) - mean[j]) * (point.get(k) - mean[k]);
                }
            }
        }
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                covarianceMatrix[j][k] /= (n - 1);
            }
        }
        return covarianceMatrix;
    }

    /**
     * Obtain the correlation coefficient matrix from the covariance matrix.
     *
     * @param covarianceMatrix The original data covariance matrix
     * @return The correlation coefficient matrix
     */
    public static double[][] calculateCorrelation(double[][] covarianceMatrix) {
        int size = covarianceMatrix.length;
        double[][] correlationMatrix = new double[size][size];
        double[] standardDeviations = new double[size];
        // 计算每个变量的标准差
        for (int i = 0; i < size; i++) {
            standardDeviations[i] = Math.sqrt(covarianceMatrix[i][i]);
        }
        // 计算相关系数矩阵
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (standardDeviations[i] == 0 || standardDeviations[j] == 0) {
                    correlationMatrix[i][j] = 0.0; // 防止除零
                } else {
                    correlationMatrix[i][j] = covarianceMatrix[i][j] /
                            (standardDeviations[i] * standardDeviations[j]);
                }
            }
        }
        return correlationMatrix;
    }

    /**
     * Find the minimum mean value (min), then scale each dimension's mean value using the formula: mean = (mean / min) * mu.
     *
     * @param mu   The mean of the Gaussian noise
     * @param mean The mean of the original data
     * @return The scaled mean after multiplication
     */
    public static double[] scaleMeanByMu(double[] mean, double mu) {
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        double[] mean2 = new double[mean.length];
        for (double num : mean) {
            if (num < min) {
                min = num;
            }
            if (num > max) {
                max = num;
            }
        }
        for (int i = 0; i < mean.length; i++) {
            mean2[i] = mu * mean[i] / min;
        }
        return mean2;
    }

    /**
     * Generate a covariance matrix for n-dimensional Gaussian noise based on the univariate variance sigmaSquare and the original data's correlation coefficient matrix.
     *
     * @param correlationMatrix The correlation coefficient matrix of the original data
     * @param sigmaSquare       The variance of the univariate Gaussian noise
     * @return The covariance matrix for the multivariate Gaussian noise
     */
    public static double[][] generateCovByCor(double[][] correlationMatrix, double sigmaSquare) {
        int size = correlationMatrix.length;
        double[][] covarianceMatrix = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                covarianceMatrix[i][j] = correlationMatrix[i][j] * sigmaSquare;
            }
        }
        return covarianceMatrix;
    }

    /**
     * 将gps的经纬度数据转化为直角坐标系XY坐标数据并缩放
     */
    public static void convertGPS(){
        Scanner scanner = new Scanner(System.in);

        System.out.print("Please input the index of the track: ");
        int trackNumber = scanner.nextInt();

        String prefix = Constants.rawPathPrefix + "gps/go_track_trackspoints_track";
        String inputCsvFile = prefix + trackNumber + ".csv"; // 输入CSV文件路径
        String outputCsvFile = prefix + trackNumber + "_raw.csv"; // 输出CSV文件路径

        List<double[]> xyCoordinates = new ArrayList<>();
        CSVFormat csvFormat = CSVFormat.DEFAULT.builder()
                .setHeader()
                .setSkipHeaderRecord(true)
                .build();
        // 读取CSV文件
        try (Reader in = new FileReader(inputCsvFile);
             CSVParser parser = new CSVParser(in, csvFormat)) {
            for (CSVRecord record : parser) {
                String la_str = record.get("latitude");
                String long_str = record.get("longitude");
                if (la_str.isEmpty()) {
                    break;
                }
                double latitude = Double.parseDouble(la_str);
                double longitude = Double.parseDouble(long_str);
                // 转换纬度和经度到x，y坐标
                double[] xy = latLongToXY(latitude, longitude);
                xyCoordinates.add(xy);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        // 对xy坐标序列进行去中心化
        double[] xMean = {0.0, 0.0};
        for (double[] xy : xyCoordinates) {
            xMean[0] += xy[0];
            xMean[1] += xy[1];
        }
        xMean[0] /= xyCoordinates.size();
        xMean[1] /= xyCoordinates.size();

        for (double[] xy : xyCoordinates) {
            xy[0] -= xMean[0];
            xy[1] -= xMean[1];
        }
        // 写入新的CSV文件
        try (Writer out = new FileWriter(outputCsvFile);
             CSVPrinter printer = new CSVPrinter(out, CSVFormat.DEFAULT)) {
            for (double[] xy : xyCoordinates) {
                printer.printRecord(xy[0], xy[1]);
            }
            System.out.println("[INFO]Done");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double[] latLongToXY(double latitude, double longitude) {
        // 米勒投影的参数和公式转换
        double R = 6378137; // 地球半径，单位米
        double L = 2 * Math.PI * R; // 地球周长
        double W = L; // 地图宽度
        double H = L / 2; // 地图高度
        double mill = 2.3; // 米勒投影参数
        double x = longitude * Math.PI / 180;
        double y = latitude * Math.PI / 180;
        y = 1.25 * Math.log(Math.tan(0.25 * Math.PI + 0.4 * y));
        // 转换到x，y坐标
        x = (W / 2) + (W / (2 * Math.PI)) * x;
        y = (H / 2) - (H / (2 * mill)) * y;
        return new double[]{x / 100, y / 100};
    }

    public static double calculateRMSE(TimeSeries ts) {
        double sumSquaredErrors = 0;
        int n = ts.size();
        int rows = ts.getDimension();
        for (int i = 0; i < n; i++) {
            for (int row = 0; row < rows; row++) {
                if(ts.getLabelAt(i)){
                    continue;
                }
                double error = ts.getTruthAt(i).get(row, 0) - ts.getRepairAt(i).get(row, 0);
                sumSquaredErrors += error * error;
            }
        }
        return Math.sqrt(sumSquaredErrors / n);
    }

    public static double calculateRMSE(UniTimeSeries ts) {
        double sumSquaredErrors = 0;
        int n = ts.size();
        for (int i = 0; i < n; i++) {
            if(ts.getLabelAt(i)){
                continue;
            }
            double error = ts.getTruthAt(i) - ts.getRepairAt(i);
            sumSquaredErrors += error * error;
        }
        return Math.sqrt(sumSquaredErrors / n);
    }

    public static double calculateRMSE_dirty(TimeSeries ts) {
        double sumSquaredErrors = 0;
        int n = ts.size();
        int rows = ts.getDimension();
        for (int i = 0; i < n; i++) {
            for (int row = 0; row < rows; row++) {
                double error = ts.getTruthAt(i).get(row, 0) - ts.getObsAt(i).get(row, 0);
                sumSquaredErrors += error * error;
            }
        }
        return Math.sqrt(sumSquaredErrors / n);
    }

    public static double calculateRMSE_dirty(UniTimeSeries ts) {
        double sumSquaredErrors = 0;
        int n = ts.size();
        for (int i = 0; i < n; i++) {
            double error = ts.getTruthAt(i) - ts.getObsAt(i);
            sumSquaredErrors += error * error;
        }
        return Math.sqrt(sumSquaredErrors / n);
    }

    public static double calculateMAE(TimeSeries ts) {
        double sumAbsoluteErrors = 0;
        int n = ts.size();
        int rows = ts.getDimension();
        for (int i = 0; i < n; i++) {
            for (int row = 0; row < rows; row++) {
                if(ts.getLabelAt(i)){
                    continue;
                }
                double error = Math.abs(ts.getTruthAt(i).get(row, 0) - ts.getRepairAt(i).get(row, 0));
                sumAbsoluteErrors += error;
            }
        }
        return sumAbsoluteErrors / n;
    }

    public static double calculateMAE(UniTimeSeries ts) {
        double sumAbsoluteErrors = 0;
        int n = ts.size();
        for (int i = 0; i < n; i++) {
            if(ts.getLabelAt(i)){
                continue;
            }
            double error = Math.abs(ts.getTruthAt(i) - ts.getRepairAt(i));
            sumAbsoluteErrors += error;
        }
        return sumAbsoluteErrors / n;
    }

    public static double calculateMAE_dirty(TimeSeries ts) {
        double sumAbsoluteErrors = 0;
        int n = ts.size();
        int rows = ts.getDimension();
        for (int i = 0; i < n; i++) {
            for (int row = 0; row < rows; row++) {
                double error = Math.abs(ts.getTruthAt(i).get(row, 0) - ts.getObsAt(i).get(row, 0));
                sumAbsoluteErrors += error;
            }
        }
        return sumAbsoluteErrors / n;
    }

    public static double calculateMAE_dirty(UniTimeSeries ts) {
        double sumAbsoluteErrors = 0;
        int n = ts.size();
        for (int i = 0; i < n; i++) {
            double error = Math.abs(ts.getTruthAt(i) - ts.getObsAt(i));
            sumAbsoluteErrors += error;
        }
        return sumAbsoluteErrors / n;
    }

    /**
     * 计算多维数据的RMSE和MAE
     *
     * @return 四元组(rmse, mae, rmse_dirty, mae_dirty)
     */
    public static List<Double> scoreResult(TimeSeries ts) {
        List<Double> score;
        // 计算RMSE和MAE作为评价指标
        double rmsDirty = calculateRMSE_dirty(ts);
        double rmse = calculateRMSE(ts);
        double maeDirty = calculateMAE_dirty(ts);
        double mae = calculateMAE(ts);
        score = Arrays.asList(rmse, mae, rmsDirty, maeDirty);
//        System.out.println("[INFO]RMSE:\t\t\t" + rmse + ".\n[INFO]MAE:\t\t\t" + mae + ".\n[INFO]RMSE_DIRTY:\t" + "" + rmsDirty + ".\n[INFO]MAE_DIRTY:\t" + maeDirty + ".\n");
        System.out.println("[INFO]RMSE:\t\t\t" + rmse + ".\n[INFO]RMSE_DIRTY:\t" + rmsDirty + ".\n");
        return score;
    }

    /**
     * 计算单维数据的RMSE和MAE
     *
     * @return 四元组(rmse, mae, rmse_dirty, mae_dirty)
     */
    public static List<Double> scoreResult(UniTimeSeries ts) {
        List<Double> score;
        double rmsDirty = calculateRMSE_dirty(ts);
        double rmse = calculateRMSE(ts);
        double maeDirty = calculateMAE_dirty(ts);
        double mae = calculateMAE(ts);
        score = Arrays.asList(rmse, mae, rmsDirty, maeDirty);
//        System.out.println("[INFO]RMSE:\t\t\t" + rmse + ".\n[INFO]MAE:\t\t\t" + mae + ".\n[INFO]RMSE_DIRTY:\t" + "" + rmsDirty + ".\n[INFO]MAE_DIRTY:\t" + maeDirty + ".\n");
        System.out.println("[INFO]RMSE:\t\t\t" + rmse + ".\n[INFO]RMSE_DIRTY:\t" + rmsDirty + ".\n");
        return score;
    }

    public static void print2DArray(double[][] array) {
        for (double[] doubles : array) {
            for (double aDouble : doubles) {
                System.out.print(aDouble + " ");
            }
            System.out.println();
        }
    }
}
