package com.MINOR.entity;

/**
 * store the result of repair
 */
public class RepairResult {
    private double rmse;
    private double timeCost;
    private int iterationNum;

    public RepairResult(double rmse, double timeCost, int iterationNum) {
        this.rmse = rmse;
        this.timeCost = timeCost;
        this.iterationNum = iterationNum;
    }

    public double getRmse() {
        return rmse;
    }

    public void setRmse(double rmse) {
        this.rmse = rmse;
    }

    public double getTimeCost() {
        return timeCost;
    }

    public void setTimeCost(double timeCost) {
        this.timeCost = timeCost;
    }

    public int getIterationNum() {
        return iterationNum;
    }

    public void setIterationNum(int iterationNum) {
        this.iterationNum = iterationNum;
    }

    /**
     * get average result
     *
     * @param seedsNum total number of seeds
     */
    public void setAve(int seedsNum) {
        this.rmse /= seedsNum;
        this.timeCost /= seedsNum;
        this.iterationNum /= seedsNum;
    }

    /**
     *  return this + other
     */
    public void add(RepairResult other) {
        this.rmse += other.rmse;
        this.timeCost += other.timeCost;
        this.iterationNum += other.iterationNum;
    }

    @Override
    public String toString() {
        return "[RESULT INFO]RepairResult{" +
                "rmse=" + rmse +
                ", timeCost=" + timeCost +
                ", iterationNum=" + iterationNum +
                '}';
    }
}
