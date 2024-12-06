package com.MINOR.Utils;

import com.MINOR.entity.TimePoint;
import com.MINOR.entity.TimeSeries;
import com.MINOR.entity.UniTimePoint;
import com.MINOR.entity.UniTimeSeries;
import org.ejml.data.DMatrixRMaj;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import static com.MINOR.Utils.MyMathUtils.scoreResult;

public class MyFileUtils {
    /**
     * Reads a CSV file into a List<List<Double>>
     *
     * @param filePath        Path to the CSV file
     * @param startColumn     Starting column
     * @param numberOfColumns Number of columns to read
     * @return Multi-variate Double sequence
     */
    public static List<List<Double>> csv2List(String filePath, int startColumn, int numberOfColumns) {
        List<List<Double>> result = new ArrayList<>();
        int length = startColumn + numberOfColumns;
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                if (values.length < startColumn + numberOfColumns) {
                    System.err.println("[ERROR]Line has insufficient columns: " + line);
                    continue;
                }
                List<Double> rowData = new ArrayList<>(numberOfColumns);
                for (int i = startColumn; i < length; i++) {
                    try {
                        rowData.add(Double.parseDouble(values[i]));
                    } catch (NumberFormatException e) {
                        System.err.println("[ERROR]Failed to parse value: " + values[i] + " at line " + result.size());
                        rowData.add(null); // or handle differently based on requirements
                    }
                }
                result.add(rowData);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return result;
    }

    /**
     * Reads a CSV file into a List<Integer>
     *
     * @param filePath    Path to the CSV file
     * @param columnIndex Starting column
     * @return Uni-variate Integer sequence
     */
    public static List<Integer> csv2IntList(String filePath, int columnIndex) {
        List<Integer> columnData = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                if (columnIndex < values.length) {
                    try {
                        columnData.add(Integer.parseInt(values[columnIndex].trim()));
                    } catch (NumberFormatException e) {
                        System.err.println("[ERROR]Invalid integer value in column " + columnIndex + ": " + values[columnIndex]);
                    }
                } else {
                    System.err.println("[WARNING]]Line has fewer columns than expected. Skipping.");
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return columnData;
    }

    /**
     * Truncates the original CSV file to the first n lines and writes to a new CSV file.
     *
     * @param inputFilePath  the path to the original CSV file
     * @param outputFilePath the path where the new CSV file will be saved
     * @param n              the number of lines to read and write
     * @throws IOException if an I/O error occurs
     */
    public static void truncateCSV(String inputFilePath, String outputFilePath, int n) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(inputFilePath));
             BufferedWriter writer = new BufferedWriter(new FileWriter(outputFilePath))) {
            String line;
            int count = 0;
            while ((line = reader.readLine()) != null && count < n) {
                writer.write(line);
                writer.newLine();
                count++;
            }
        }
    }

    public static void csvToListInt(String csvFilePath, List<Integer> t, int col) {
        try (BufferedReader br = new BufferedReader(new FileReader(csvFilePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                if (values[col].matches("-?\\d+(\\.\\d+)?([eE][-+]?\\d+)?")) {
                    t.add(Integer.parseInt(values[col]));
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Reads two CSV files and converts them into a UniTimeSeries object
     *
     * @param ts        Time series data object
     * @param dataFPath Path to the data CSV file
     */
    public static void csv2UTS(UniTimeSeries ts, String dataFPath) {
        try (BufferedReader br = new BufferedReader(new FileReader(dataFPath))
        ) {
            ts.s.clear();
            ts.clearCandidates();
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                if (values.length >= 3) { // +1 for timestamp
                    long timestamp = Long.parseLong(values[0]); // timestamp
                    // initialize observe data and ground truth
                    double obsData, truthData;
                    boolean lb;
                    // Parse the 2nd column as valObs and the 3rd column as valTruth.
                    if (values[1].matches("-?\\d+(\\.\\d+)?([eE][-+]?\\d+)?") &&
                            values[2].matches("-?\\d+(\\.\\d+)?([eE][-+]?\\d+)?")) {
                        obsData = Double.parseDouble(values[1]);
                        truthData = Double.parseDouble(values[2]);
                        lb = MyMathUtils.int2Bool(Integer.parseInt(values[3]));
                    } else {
                        break;
                    }
                    double valObs = obsData;
                    double valTruth = truthData;
                    UniTimePoint timePoint = new UniTimePoint(timestamp, valTruth, valObs, Double.MAX_VALUE, valObs, lb);
                    ts.s.add(timePoint);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Reads two CSV files and converts them into a TimeSeries object with the corresponding dimension
     *
     * @param ts        Time series data object
     * @param dimension Dimension of the time series data
     * @param dataFPath Path to the data CSV file
     */
    public static void csv2MTS(TimeSeries ts, int dimension, String dataFPath) {
        try (BufferedReader br = new BufferedReader(new FileReader(dataFPath))) {
            ts.s.clear();
            ts.clearCandidates();
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                if (values.length >= dimension * 2 + 2) { // +2 for timestamp and label
                    long timestamp = Long.parseLong(values[0]); // timestamp
                    double[] obsData = new double[dimension];
                    double[] truthData = new double[dimension];
                    boolean lb;
                    boolean isValid = true;
                    // Parse the first m columns as valObs and the last m columns as valTruth.
                    for (int i = 0; i < dimension; i++) {
                        if (values[i + 1].matches("-?\\d+(\\.\\d+)?([eE][-+]?\\d+)?") &&
                                values[i + dimension + 1].matches("-?\\d+(\\.\\d+)?([eE][-+]?\\d+)?")) {
                            obsData[i] = Double.parseDouble(values[i + 1]);
                            truthData[i] = Double.parseDouble(values[i + dimension + 1]);
                        } else {
                            isValid = false;
                            break;
                        }
                    }
                    lb = MyMathUtils.int2Bool(Integer.parseInt(values[dimension * 2 + 1]));
                    if (isValid) {
                        // The generated output is a column vector.
                        DMatrixRMaj valObs = new DMatrixRMaj(obsData);
                        DMatrixRMaj valTruth = new DMatrixRMaj(truthData);
                        DMatrixRMaj valRepaired = new DMatrixRMaj(obsData);
                        TimePoint timePoint = new TimePoint(timestamp, dimension, valTruth, valObs, null, valRepaired, lb);
                        ts.s.add(timePoint);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Calculate RMSE and other metrics for the repaired univariate time series data, and save the repaired data to a CSV file.
     *
     * @param csvFilePath The file path; if null, the data will not be output to a file, only the scores will be output.
     * @param ts          The repaired time series data.
     */
    public static void UTSExporter(String csvFilePath, UniTimeSeries ts) {
        List<Double> score = scoreResult(ts);
        if (csvFilePath == null) {
            System.out.println("[INFO]null resultCsvFilePath.");
            return;
        }
        try (FileWriter writer = new FileWriter(csvFilePath)) {
            int n = ts.size();
            assert (n >= 4);
            for (int i = 0; i < n; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(ts.getTimestampAt(i));
                sb.append(",").append(ts.getObsAt(i));
                sb.append(",").append(ts.getRepairAt(i));
                sb.append(",").append(ts.getTruthAt(i));
                sb.append(",").append(MyMathUtils.bool2Int(ts.getLabelAt(i)));
                switch (i) {
                    case 0:
                        sb.append(",RMSE,").append(score.get(i));
                        break;
                    case 1:
                        sb.append(",MAE,").append(score.get(i));
                        break;
                    case 2:
                        sb.append(",RMSE_Dirty,").append(score.get(i));
                        break;
                    case 3:
                        sb.append(",MAE_Dirty,").append(score.get(i));
                        break;
                    default:
                }
                sb.append('\n');
                writer.write(sb.toString());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Calculate RMSE and other metrics for the repaired multivariate time series data, and save the repaired data to a CSV file.
     *
     * @param csvFilePath The file path; if null, the data will not be output to a file, only the scores will be output.
     * @param ts          The repaired time series data.
     */
    public static void MTSExporter(String csvFilePath, TimeSeries ts) {
        List<Double> score = scoreResult(ts);
        if (csvFilePath == null) {
            System.out.println("[INFO]null resultCsvFilePath.");
            return;
        }
        try (FileWriter writer = new FileWriter(csvFilePath)) {
            int n = ts.size();
            int m = ts.getDimension();
            assert (n >= 4);
            for (int i = 0; i < n; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(ts.getTimestampAt(i));
                for (int j = 0; j < m; j++) {
                    sb.append(",").append(ts.getObsAt(i).get(j, 0));
                }
                for (int j = 0; j < m; j++) {
                    sb.append(",").append(ts.getRepairAt(i).get(j, 0));
                }
                for (int j = 0; j < m; j++) {
                    sb.append(",").append(ts.getTruthAt(i).get(j, 0));
                }
                sb.append(",").append(MyMathUtils.bool2Int(ts.getLabelAt(i)));
                switch (i) {
                    case 0:
                        sb.append(",RMSE,").append(score.get(i));
                        break;
                    case 1:
                        sb.append(",MAE,").append(score.get(i));
                        break;
                    case 2:
                        sb.append(",RMSE_Dirty,").append(score.get(i));
                        break;
                    case 3:
                        sb.append(",MAE_Dirty,").append(score.get(i));
                        break;
                    default:
                }
                sb.append('\n');
                writer.write(sb.toString());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Write the phi calculated by the model to a file, applicable for both univariate and multivariate model phi.
     *
     * @param csvFilePath The path to the file where phi will be stored.
     * @param phi         The phi parameter of the model.
     */
    public static void PhiExporter(String csvFilePath, DMatrixRMaj phi) {
        if (csvFilePath == null) {
            return;
        }
        try (FileWriter writer = new FileWriter(csvFilePath)) {
            int nr = phi.numRows;
            int nc = phi.numCols;
            for (int i = 0; i < nr; i++) {
                StringBuilder sb = new StringBuilder();
                for (int j = 0; j < nc; j++) {
                    sb.append(phi.get(i, j));
                    if (j < nc - 1) {
                        sb.append(",");
                    }
                }
                sb.append('\n');
                writer.write(sb.toString());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * export data of application data set
     *
     * @param outputFile    file name to save data
     * @param timestamps    timestamps
     * @param datas         dirty or reapaired data
     * @throws IOException  -
     */
    public static void exportAppData(String outputFile, List<Long> timestamps, List<List<Double>> datas) throws IOException {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile))) {
            int numLines = timestamps.size();
            int numColumns = datas.size();

            for (int i = 0; i < numLines; i++) {
                StringBuilder line = new StringBuilder();
                line.append(timestamps.get(i));
                for (int j = 0; j < numColumns; j++) {
                    line.append(",").append(datas.get(j).get(i));
                }
                bw.write(line.toString());
                bw.newLine();
            }
        }
    }

    /**
     * export labels for MINOR-B, MINOR-U, VARX, IMR to repair data of application data set
     *
     * @param labelsFile    file name to save label
     * @param dataSize      size of data
     * @param labels        a set of labels
     * @throws IOException  -
     */
    public static void exportAppLabels(String labelsFile, int dataSize, List<Set<Integer>> labels) throws IOException {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(labelsFile))) {
            int dim = labels.size();
            for (int i = 0; i < dataSize; i++) {
                StringBuilder sb = new StringBuilder();
                for (int j = 0; j < dim; j++) {
                    sb.append(labels.get(j).contains(i) ? 1 : 0);
                    sb.append(",");
                }
                if (!sb.isEmpty()) {
                    sb.deleteCharAt(sb.length() - 1);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
        }
    }

    public static List<List<Double>> readAppData(String csvFile){
        List<List<Double>> datas = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                // Remaining columns as data
                for (int i = 1; i < values.length; i++) {
                    if (datas.size() < i) {
                        datas.add(new ArrayList<>());
                    }
                    datas.get(i - 1).add(Double.parseDouble(values[i]));
                }
            }
        }catch (IOException e) {
            e.printStackTrace();
        }
        return datas;
    }

    public static List<List<Integer>> readAppLabel(String labelFile) {
        List<List<Integer>> columnData = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(labelFile))) {
            String line;
            List<String[]> allRows = new ArrayList<>();
            // 读取每一行并解析其为整数数组
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                allRows.add(values);
            }
            if (allRows.isEmpty()) {
                return columnData; // 如果文件为空或没有有效的数据，返回空的列表
            }
            // 初始化列
            int numColumns = allRows.get(0).length;
            for (int i = 0; i < numColumns; i++) {
                columnData.add(new ArrayList<>());
            }
            // 填充列数据
            for (String[] row : allRows) {
                for (int i = 0; i < numColumns; i++) {
                    try {
                        columnData.get(i).add(Integer.parseInt(row[i].trim()));
                    } catch (NumberFormatException e) {
                        // 处理无法解析为整数的情况
                        System.err.println("[ERROR]invalid label: " + row[i]);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return columnData;
    }

    /**
     * write RMSE or time cost to csv files
     *
     * @param fileName   csv file name
     * @param modelNames names of  models used in experiment
     * @param varVals    variable values of the experiment
     * @param costs      rmse or time cost
     */
    public static void exportResult(String fileName, String[] modelNames, double[] varVals, double[][] costs) {
        try {
            FileWriter fw = new FileWriter(fileName);
            PrintWriter pw = new PrintWriter(fw);
            int SI = 0;
            if (modelNames[0].equals(" ")) {
                SI = 1;
            }
            int size = varVals.length;
            int modelCnt = costs[0].length;
            pw.print(" ");
            for (double varVal : varVals) {
                pw.print(",");
                pw.print(varVal);
            }
            pw.println();
            for (int j = 0; j < modelCnt; j++) {
                pw.print(modelNames[j + SI]);
                for (int i = 0; i < size; i++) {
                    pw.print(",");
                    pw.print(costs[i][j]);
                }
                pw.println();
            }
            pw.flush();
            pw.close();
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void exportResult(String fileName, String[] modelNames, int[] varVals, double[][] costs) {
        try {
            FileWriter fw = new FileWriter(fileName);
            PrintWriter pw = new PrintWriter(fw);
            int SI = 0;
            if (modelNames[0].equals(" ")) {
                SI = 1;
            }
            int size = varVals.length;
            int modelCnt = costs[0].length;
            pw.print(" ");
            for (int varVal : varVals) {
                pw.print(",");
                pw.print(varVal);
            }
            pw.println();
            for (int j = 0; j < modelCnt; j++) {
                pw.print(modelNames[j + SI]);
                for (int i = 0; i < size; i++) {
                    pw.print(",");
                    pw.print(costs[i][j]);
                }
                pw.println();
            }
            pw.flush();
            pw.close();
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void exportResult(String fileName, String[] modelNames, double[] varVals, int[][] costs) {
        try {
            FileWriter fw = new FileWriter(fileName);
            PrintWriter pw = new PrintWriter(fw);
            int SI = 0;
            if (modelNames[0].equals(" ")) {
                SI = 1;
            }
            int size = varVals.length;
            int modelCnt = costs[0].length;
            pw.print(" ");
            for (double varVal : varVals) {
                pw.print(",");
                pw.print(varVal);
            }
            pw.println();
            for (int j = 0; j < modelCnt; j++) {
                pw.print(modelNames[j + SI]);
                for (int i = 0; i < size; i++) {
                    pw.print(",");
                    pw.print(costs[i][j]);
                }
                pw.println();
            }
            pw.flush();
            pw.close();
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void exportResult(String fileName, String[] modelNames, int[] varVals, int[][] costs) {
        try {
            FileWriter fw = new FileWriter(fileName);
            PrintWriter pw = new PrintWriter(fw);
            int SI = 0;
            if (modelNames[0].equals(" ")) {
                SI = 1;
            }
            int size = varVals.length;
            int modelCnt = costs[0].length;
            pw.print(" ");
            for (int varVal : varVals) {
                pw.print(",");
                pw.print(varVal);
            }
            pw.println();
            for (int j = 0; j < modelCnt; j++) {
                pw.print(modelNames[j + SI]);
                for (int i = 0; i < size; i++) {
                    pw.print(",");
                    pw.print(costs[i][j]);
                }
                pw.println();
            }
            pw.flush();
            pw.close();
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
