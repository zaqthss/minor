package com.MINOR.Utils;

import com.MINOR.entity.UniTimePoint;
import com.MINOR.entity.UniTimeSeries;

import java.util.List;

/**
 * Utils for Akane Heuristic
 */
public class AkaneUtils {
    public static class Triplet implements Comparable<Triplet> {
        private int first;
        private String second;
        private Double third;

        public Triplet(int first, String second, Double third) {
            this.first = first;
            this.second = second;
            this.third = third;
        }

        public int getFirst() {
            return first;
        }

        public void setFirst(int first) {
            this.first = first;
        }

        public String getSecond() {
            return second;
        }

        public void setSecond(String second) {
            this.second = second;
        }

        public Double getThird() {
            return third;
        }

        public void setThird(Double third) {
            this.third = third;
        }

        @Override
        public int compareTo(Triplet other) {
            return this.third.compareTo(other.getThird());
        }

        @Override
        public String toString() {
            return "(" + first + ", " + second + ", " + third + ")";
        }
    }

    public void heapify(Triplet arr[], int N, int i, int[] posRecord) {
        int largest = i; // Initialize largest as root
        int l = 2 * i + 1; // left = 2*i + 1
        int r = 2 * i + 2; // right = 2*i + 2

        // If left child is larger than root
        if (l < N && arr[l].compareTo(arr[largest]) > 0)
            largest = l;

        // If right child is larger than largest so far
        if (r < N && arr[r].compareTo(arr[largest]) > 0)
            largest = r;

        // If largest is not root
        if (largest != i) {
            Triplet swap = arr[i];
            arr[i] = arr[largest];
            arr[largest] = swap;

            int tempPos = posRecord[arr[i].getFirst()];
            posRecord[arr[i].getFirst()] = posRecord[arr[largest].getFirst()];
            posRecord[arr[largest].getFirst()] = tempPos;

            // Recursively heapify the affected sub-tree
            heapify(arr, N, largest, posRecord);
        }
    }

    // Function to build a Max-Heap from the Array
    public void buildHeap(Triplet[] arr, int[] posRecord, int N) {
        // Index of last non-leaf node
        int startIdx = (N / 2) - 1;

        // Perform reverse level order traversal
        // from last non-leaf node and heapify
        // each node
        for (int i = startIdx; i >= 0; i--) {
            heapify(arr, N, i, posRecord);
        }
    }

    public void replaceHeap(Triplet[] arr, int[] posRecord, int index, Triplet sub, int N) {
        int targetIndex = posRecord[index];
        // Step 2: Replace the target element with the given Triplet sub
        arr[targetIndex] = sub;

        // Step 3: Reheapify the heap after replacement
        // Check and adjust upwards
        while (targetIndex != 0 && arr[targetIndex].compareTo(arr[(targetIndex - 1) / 2]) > 0) {
            // Swap current node with its parent
            Triplet temp = arr[targetIndex];
            arr[targetIndex] = arr[(targetIndex - 1) / 2];
            arr[(targetIndex - 1) / 2] = temp;

            int tempPos = posRecord[arr[targetIndex].getFirst()];
            posRecord[arr[targetIndex].getFirst()] = posRecord[arr[(targetIndex - 1) / 2].getFirst()];
            posRecord[arr[(targetIndex - 1) / 2].getFirst()] = tempPos;

            // Move to the parent index
            targetIndex = (targetIndex - 1) / 2;
        }

        // Check and adjust downwards
        heapify(arr, N, targetIndex, posRecord);
    }

    // A utility function to print the array
    // representation of Heap
    public void printHeap(Triplet arr[], int N) {
        System.out.println(
                "Array representation of Heap is:");

        for (int i = 0; i < N; ++i)
            System.out.println(arr[i].toString());

        System.out.println();
    }

    public double[] TS2List(UniTimeSeries timeseries) {
        List<UniTimePoint> totalList = timeseries.s;
        int n = totalList.size();

        double[] values = new double[n];

        for (int i = 0; i < n; i++) {
            values[i] = totalList.get(i).getValObs();
        }

        return values;
    }

    public UniTimeSeries List2TS(double[] values) {
        int n = values.length;

        UniTimeSeries resultSeries = new UniTimeSeries();
        for (int i = 0; i < n; i++) {
            UniTimePoint tp = new UniTimePoint(i + 1, values[i]);
            resultSeries.add(tp);
        }

        return resultSeries;
    }

    /**
     * Apply the repair results to the original data.
     *
     * @param clean  cleaned data without truth and labels
     * @param origin original input time series
     * @return time series with repaired data
     */
    public UniTimeSeries applyCleanedData(double[] clean, UniTimeSeries origin) {
        int n = clean.length;
        for (int i = 0; i < n; i++) {
            origin.setRepairAt(i, clean[i]);
        }
        return origin;
    }
}
