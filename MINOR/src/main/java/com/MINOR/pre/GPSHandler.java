package com.MINOR.pre;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.csv.CSVRecord;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * 预处理GPS数据
 * 根据输入的标记率，去除raw data中部分标记，形成dirty data
 */
public class GPSHandler {
    private final String inputFPath;
    private final String outputFPath;
    private final double lr;  // label rate
    private final long seed;

    public GPSHandler(String inputFPath, String outputFPath, double lr, long seed) {
        this.inputFPath = inputFPath;
        this.outputFPath = outputFPath;
        this.lr = lr;
        this.seed = seed;
    }

    public void run(){
        int labelCol = 5;   // label在第6列
        try {
            // 读取CSV文件
            Reader reader = Files.newBufferedReader(Paths.get(inputFPath));
            CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT);

            List<CSVRecord> records = csvParser.getRecords();
            List<String[]> modifiedRecords = new ArrayList<>();

            int numCol = records.size();
            int numLabels = 0;
            List<Integer> lbRows = new ArrayList<>();
            List<Integer> keptRows;

            for (int i = 0; i < numCol; i++) {
                if (records.get(i).get(labelCol).equals("1")) { // 检查第6列是否为1
                    numLabels++;
                    lbRows.add(i);
                }
            }
            int numToKeep = (int) Math.ceil(numLabels * lr);
            Random rand = new Random(seed);
            Collections.shuffle(lbRows, rand);
            keptRows = lbRows.subList(0, numToKeep);    // [0,numToKeep-1)

            for (int i = 0; i < numCol; i++) {
                String[] modifiedRecord = new String[7];
                CSVRecord record = records.get(i);
                // 复制前5列
                for (int j = 0; j < 5; j++) {
                    modifiedRecord[j] = record.get(j);
                }
                // 处理第6列
                if (record.get(labelCol).equals("1") && keptRows.contains(i)) {
                    modifiedRecord[5] = "1";
                } else {
                    modifiedRecord[5] = "0";
                }

                modifiedRecords.add(modifiedRecord);
            }
            csvParser.close();
            // 写入修改后的CSV文件
            BufferedWriter writer = Files.newBufferedWriter(Paths.get(outputFPath));
            CSVPrinter csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT);

            for (String[] modifiedRecord : modifiedRecords) {
                csvPrinter.printRecord((Object[]) modifiedRecord);
            }
            csvPrinter.flush();
            csvPrinter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
