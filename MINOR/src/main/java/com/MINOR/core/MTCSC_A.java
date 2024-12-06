package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.RepairResult;
import com.MINOR.entity.TimePoint;
import com.MINOR.entity.TimeSeries;
import org.ejml.data.DMatrixRMaj;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Multivariate Time Series Cleaning under Speed Constraints
 * MTCSC-A
 */
public class MTCSC_A {

    private TimeSeries timeseries;
    private TimePoint kp;
    private int T;       // the window size
    private double SMAX;  // maximum speed
    // private ArrayList<Double> SpeedList = new ArrayList<Double>(); // speed windows
    // private ArrayList<Double> SortedSpeedList = new ArrayList<Double>(); // speed windows
    private ArrayList<Double> preSpeedList = new ArrayList<Double>();
    private ArrayList<Double> nowSpeedList = new ArrayList<Double>();
    private ArrayList<Integer> preSpeedDistribution = new ArrayList<Integer>();
    private ArrayList<Integer> nowSpeedDistribution = new ArrayList<Integer>();
    // private int Index = 95; // speed rate
    private double Drate = 0.025;
    private double Beta = 0.95;
    private double Threshold = 1.0;
    private int SWSize = 200;
    private int Bucket = 5;

    private int dimension = 2;

    /**
     * @param timeseries timeseries
     * @param sMax       maximum allowed speed
     * @param t          the window size
     */
    public MTCSC_A(TimeSeries timeseries, double sMax, int t, double drate, double threshold, int swsize, double beta, int bucket) {
        setTimeSeries(timeseries);
        setSMAX(sMax);
        setT(t);
        // setIndex(drate);
        setDrate(drate);
        setBeta(beta);
        setThreshold(threshold);
        setSWSize(swsize);
        setBucket(bucket);
        setSpeedList();
    }

    public void setTimeSeries(TimeSeries timeSeries) {
        this.timeseries = timeSeries;
    }

    public void setT(int t) {
        this.T = t;
    }

    public void setSMAX(double SMAX) {
        this.SMAX = SMAX;
    }

    public void setDrate(double drate) {
        this.Drate = drate;
    }

    public void setBeta(double beta) {
        this.Beta = beta;
    }

    public void setBucket(int bucket) {
        this.Bucket = bucket;
    }

    public void setThreshold(double threshold) {
        this.Threshold = threshold;
    }

    public void setSWSize(int swsize) {
        this.SWSize = swsize;
    }

    public void setSpeedList() {
        for (int i = 0; i < (this.Bucket + 1); i++) {
            this.preSpeedDistribution.add(0);
            this.nowSpeedDistribution.add(0);
        }
    }

    // Update speed based on KL divergence
    public void addSpeed(TimePoint prePoint, TimePoint kpPoint) {
//        double[] observe1 = prePoint.getOrgval();
//        double[] observe2 = kpPoint.getOrgval();
        DMatrixRMaj observe1 = prePoint.getValObs();
        DMatrixRMaj observe2 = kpPoint.getValObs();
        double tempSpeed = Math.sqrt(Math.pow(observe1.get(0, 0) - observe2.get(0, 0), 2) + Math.pow(observe1.get(1, 0) - observe2.get(1, 0), 2)) / (kpPoint.getTimestamp() - prePoint.getTimestamp());
        int distribution = 0;
        // add point to preSpeedList and update the preSpeedDistribution
        if (this.preSpeedList.size() < this.SWSize) {
            this.preSpeedList.add(tempSpeed);
            distribution = (int) Math.ceil(tempSpeed / this.SMAX * this.Bucket);
            if (0 < distribution && distribution <= this.Bucket) {
                int temp = this.preSpeedDistribution.get(distribution - 1) + 1;
                this.preSpeedDistribution.set(distribution - 1, temp);
            } else if (distribution == 0) {
                int temp = this.preSpeedDistribution.get(0) + 1;
                this.preSpeedDistribution.set(0, temp);
            } else {
                int temp = this.preSpeedDistribution.get(this.Bucket) + 1;
                this.preSpeedDistribution.set(this.Bucket, temp);
            }
        }
        // add point to nowSpeedList and update the nowSpeedDistribution
        else if (this.preSpeedList.size() == this.SWSize && this.nowSpeedList.size() < this.SWSize) {
            this.nowSpeedList.add(tempSpeed);
            distribution = (int) Math.ceil(tempSpeed / this.SMAX * this.Bucket);
            if (0 < distribution && distribution <= this.Bucket) {
                int temp = this.nowSpeedDistribution.get(distribution - 1) + 1;
                this.nowSpeedDistribution.set(distribution - 1, temp);
            } else if (distribution == 0) {
                int temp = this.nowSpeedDistribution.get(0) + 1;
                this.nowSpeedDistribution.set(0, temp);
            } else {
                int temp = this.nowSpeedDistribution.get(this.Bucket) + 1;
                this.nowSpeedDistribution.set(this.Bucket, temp);
            }
        } else if (this.preSpeedList.size() == this.SWSize && this.nowSpeedList.size() == this.SWSize) {
            // calculate the KL between preSpeedDistribution and nowSpeedDistribution
            double klDiv = KL();
            // remove the first point of the preSpeedList
            double firstSpeed = this.preSpeedList.get(0);
            distribution = (int) Math.ceil(firstSpeed / this.SMAX * this.Bucket);
            if (0 < distribution && distribution <= this.Bucket) {
                int temp = this.preSpeedDistribution.get(distribution - 1) - 1;
                this.preSpeedDistribution.set(distribution - 1, temp);
            } else if (distribution == 0) {
                int temp = this.preSpeedDistribution.get(0) - 1;
                this.preSpeedDistribution.set(0, temp);
            } else {
                int temp = this.preSpeedDistribution.get(this.Bucket) - 1;
                this.preSpeedDistribution.set(this.Bucket, temp);
            }
            this.preSpeedList.remove(0);
            firstSpeed = this.nowSpeedList.get(0);
            this.preSpeedList.add(firstSpeed);
            distribution = (int) Math.ceil(firstSpeed / this.SMAX * this.Bucket);
            if (0 < distribution && distribution <= this.Bucket) {
                // add a point to the preSpeedList
                int temp = this.preSpeedDistribution.get(distribution - 1) + 1;
                this.preSpeedDistribution.set(distribution - 1, temp);
                // remove the first point of the nowSpeedList
                temp = this.nowSpeedDistribution.get(distribution - 1) - 1;
                this.nowSpeedDistribution.set(distribution - 1, temp);
            } else if (distribution == 0) {
                int temp = this.preSpeedDistribution.get(0) + 1;
                this.preSpeedDistribution.set(0, temp);
                // remove the first point of the nowSpeedList
                temp = this.nowSpeedDistribution.get(0) - 1;
                this.nowSpeedDistribution.set(0, temp);
            } else {
                int temp = this.preSpeedDistribution.get(this.Bucket) + 1;
                this.preSpeedDistribution.set(this.Bucket, temp);
                temp = this.nowSpeedDistribution.get(this.Bucket) - 1;
                this.nowSpeedDistribution.set(this.Bucket, temp);
            }
            this.nowSpeedList.remove(0);
            // add a point to the nowSpeedList
            this.nowSpeedList.add(tempSpeed);
            distribution = (int) Math.ceil(tempSpeed / this.SMAX * this.Bucket);
            if (0 < distribution && distribution <= this.Bucket) {
                int temp = this.nowSpeedDistribution.get(distribution - 1) + 1;
                this.nowSpeedDistribution.set(distribution - 1, temp);
            } else if (distribution == 0) {
                int temp = this.nowSpeedDistribution.get(0) + 1;
                this.nowSpeedDistribution.set(0, temp);
            } else {
                int temp = this.nowSpeedDistribution.get(this.Bucket) + 1;
                this.nowSpeedDistribution.set(this.Bucket, temp);
            }
            // beyond the threshold so update the speed
            if (klDiv > this.Threshold) {
                updateSpeed();
                upSpeedDistribution();
            }
            // System.out.println("KL is " + klDiv);
        } else {
            System.out.println("error!!");
        }
    }

    public double KL() {
        // calculate frequency of preSpeedDistribution
        ArrayList<Double> preFrequency = new ArrayList<>();
        for (int count : preSpeedDistribution) {
            double frequency = (double) count / (double) this.SWSize; // calculate the frequence
            preFrequency.add(frequency); // add frequence
        }

        // calculate frequency of nowSpeedDistribution
        ArrayList<Double> nowFrequency = new ArrayList<>();
        for (int count : nowSpeedDistribution) {
            double frequency = (double) count / (double) this.SWSize; // calculate the frequence
            nowFrequency.add(frequency); // add frequence
        }

        // calculate the KL between preSpeedDistribution and nowSpeedDistribution
        double klDiv = 0.0;
        double epsilon = 1 / (double) this.SWSize; // if q==0
        for (int i = 0; i < (this.Bucket + 1); i++) {
            double pVal = preFrequency.get(i);
            double qVal = nowFrequency.get(i);
            if (pVal > 0) {
                if (qVal != 0) {
                    klDiv += pVal * Math.log(pVal / qVal);
                } else {
                    klDiv += pVal * Math.log(pVal / epsilon);
                }
            }
        }
        return klDiv;
    }

    public void updateSpeed() {
        ArrayList<Double> tempSpeedList = new ArrayList<>(this.nowSpeedList);
        Collections.sort(tempSpeedList);
        // double rate = 1-2*this.Drate-0.05;
        double rate = 1 - 2 * this.Drate;
        int len = tempSpeedList.size();
        int index = (int) (len * rate);
        // double s = tempSpeedList.get(index-1) *1.3;
        double s = tempSpeedList.get(index - 1);
        s = s / this.Beta;
        this.SMAX = s;
    }

    public void upSpeedDistribution() {
        for (int i = 0; i < (this.Bucket + 1); i++) {
            this.preSpeedDistribution.set(i, 0);
            this.nowSpeedDistribution.set(i, 0);
        }
        for (int i = 0; i < this.preSpeedList.size(); i++) {
            // update the preSpeedDistribution
            double tempSpeed = this.preSpeedList.get(i);
            int distribution = (int) Math.ceil(tempSpeed / this.SMAX * this.Bucket);
            if (0 < distribution && distribution <= this.Bucket) {
                int temp = this.preSpeedDistribution.get(distribution - 1) + 1;
                this.preSpeedDistribution.set(distribution - 1, temp);
            } else {
                int temp = this.preSpeedDistribution.get(this.Bucket) + 1;
                this.preSpeedDistribution.set(this.Bucket, temp);
            }
            // update the nowSpeedDistribution
            tempSpeed = this.nowSpeedList.get(i);
            distribution = (int) Math.ceil(tempSpeed / this.SMAX * this.Bucket);
            if (0 < distribution && distribution <= this.Bucket) {
                int temp = this.nowSpeedDistribution.get(distribution - 1) + 1;
                this.nowSpeedDistribution.set(distribution - 1, temp);
            } else {
                int temp = this.nowSpeedDistribution.get(this.Bucket) + 1;
                this.nowSpeedDistribution.set(this.Bucket, temp);
            }
        }
    }

    /**
     * @return timeseries after repair
     */
    public RepairResult mainScreen() {
        long startTime = System.nanoTime();
        List<TimePoint> totalList = timeseries.s;
        int size = totalList.size();

        long preEnd = -1, curEnd;
        // the startTime in the window, the real end time in the window, the maximum allowed
        long wStartTime, wEndTime, wGoalTime;
        long curTime = 0;
        TimePoint prePoint = null;    // the last fixed point
        TimePoint tp;

        TimeSeries tempSeries = new TimeSeries();
        List<TimePoint> tempList;

        int readIndex = 1; // the point should be read in

        // initial
        tp = totalList.get(0);
        tempSeries.add(tp);
        wStartTime = tp.getTimestamp();
        wEndTime = wStartTime;
        wGoalTime = wStartTime + T;

        while (readIndex < size) {
            tp = totalList.get(readIndex);
            curTime = tp.getTimestamp();

            // This point shouldn't be added until the repair is over
            if (curTime > wGoalTime) {
                while (true) {
                    tempList = tempSeries.s;
                    if (tempList.size() == 0) {
                        // if all the points in tempList has been handled
                        tempSeries.add(tp);  // the current point should be a new start
                        // prePoint = tp;
                        wGoalTime = curTime + T;
                        wEndTime = curTime;
                        break;
                    }

                    kp = tempList.get(0);
                    wStartTime = kp.getTimestamp();
                    wGoalTime = wStartTime + T;

                    if (curTime <= wGoalTime) {
                        // then should read in new points
                        TimePoint tmpPoint = tempSeries.s.get(tempSeries.s.size() - 1);
                        addSpeed(tmpPoint, tp);
                        tempSeries.add(tp);
                        wEndTime = curTime;
                        break;
                    }

                    curEnd = wEndTime;

                    if (preEnd == -1) {
                        prePoint = kp;

                    }
                    local(tempSeries, prePoint);
                    prePoint = kp;
                    preEnd = curEnd;

                    // remove the keyPoint
                    tempSeries.s.remove(0);
                } // end of while(true)
            } else {
                if (curTime > wEndTime) {
                    // suppose the sequence is in order, so it must happen
                    TimePoint tmpPoint = tempSeries.s.get(tempSeries.s.size() - 1);
                    addSpeed(tmpPoint, tp);
                    tempSeries.add(tp);
                    wEndTime = curTime;
                }
            }

            readIndex++;  // read another one
        }

        // handle the last window
        while (tempSeries.size() > 0) {
            tempList = tempSeries.s;
            kp = tempList.get(0);
            if (prePoint == null) {
                prePoint = kp;
            }
            local(tempSeries, prePoint);
            prePoint = kp;
            tempList.remove(0);
        }

//        return timeseries;
        long endTime = System.nanoTime();
        double duration = ((double) (endTime - startTime)) / 1_000_000.0d;
        return new RepairResult(MyMathUtils.calculateRMSE(timeseries), duration, 1);
    }

    // calculate the distance
    private double distance(TimePoint prePoint, TimePoint kp) {
        double distance = 0;
//        double[] xy_pp = new double[2];
//        double[] xy_kp = new double[2];
//        xy_pp = prePoint.getModify();
//        xy_kp = kp.getModify();
        DMatrixRMaj xy_pp = prePoint.getValRepaired();
        DMatrixRMaj xy_kp = kp.getValRepaired();
        distance = (xy_pp.get(0,0) - xy_kp.get(0,0)) * (xy_pp.get(0,0)- xy_kp.get(0,0))
                + (xy_pp.get(1,0) - xy_kp.get(1,0)) * (xy_pp.get(1,0) - xy_kp.get(1,0));
        return distance;
    }

    public static double calculateDistance(TimePoint tp1, TimePoint tp2) {
//        double[] xy1 = tp1.getModify();
//        double[] xy2 = tp2.getModify();
        DMatrixRMaj xy1 = tp1.getValRepaired();
        DMatrixRMaj xy2 = tp2.getValRepaired();
        return Math.sqrt(Math.pow(xy1.get(0,0) - xy2.get(0,0), 2) + Math.pow(xy1.get(1,0) - xy2.get(1,0), 2));
    }

    // judge whether to modify
    private boolean judgeModify(DMatrixRMaj preVal, DMatrixRMaj maxVal, DMatrixRMaj kpVal,
                                long preTime, long maxTime, long kpTime) {
        boolean judge = true;
        double c1, c2;
        c1 = (preVal.get(0,0) - kpVal.get(0,0)) * (preVal.get(0,0) - kpVal.get(0,0))
                + (preVal.get(1,0) - kpVal.get(1,0)) * (preVal.get(1,0) - kpVal.get(1,0));
        c2 = (maxVal.get(0,0) - kpVal.get(0,0)) * (maxVal.get(0,0) - kpVal.get(0,0))
                + (maxVal.get(1,0) - kpVal.get(1,0)) * (maxVal.get(1,0) - kpVal.get(1,0));
        if (c1 <= (kpTime - preTime) * (kpTime - preTime) * SMAX * SMAX &&
                c2 <= (maxTime - kpTime) * (maxTime - kpTime) * SMAX * SMAX) {
            judge = false;
        }
        return judge;
    }

    /**
     * @param timeSeries timeseries in a window
     * @param prePoint   the former modified point
     */
    private void local(TimeSeries timeSeries, TimePoint prePoint) {
        List<TimePoint> tempList = timeSeries.s;

        // get bound
        long preTime = prePoint.getTimestamp();
//        double[] preVal = new double[2];
//        preVal = prePoint.getModify();
        long kpTime = kp.getTimestamp();
//        double[] kpVal = new double[2];
//        kpVal = kp.getModify();
        DMatrixRMaj preVal = prePoint.getValRepaired();
        DMatrixRMaj kpVal = kp.getValRepaired();

        // form candidates
        int length = tempList.size();
        int[] top = new int[length + 1]; //Default is 0,-1
        int[] len = new int[length + 1]; //Default is 0, take the longest to form candidate

        TimePoint tp1, tp2;
        if (length == 1) {
            if (judgeModify(preVal, preVal, kpVal, preTime, preTime, kpTime)) {
//                kp.setModify(preVal);
                kp.setValRepaired(new DMatrixRMaj(preVal));
            }
            return;
        }
        // if(kpTime == 937){
        //   double a = 0;
        // }
        int topIndex = 0;
        // find the first top
        for (int i = 1; i < length; ++i) {
            tp1 = tempList.get(i);
            long t1 = tp1.getTimestamp();
            if (distance(prePoint, tp1) <= ((t1 - preTime) * (t1 - preTime) * SMAX * SMAX)) {
                top[i] = -1;
                len[i] = 1;
                topIndex = i;
                break;
            }
        }
        // handle xk+topIndex+1 to xk+w
        boolean flag = false;
        for (int i = topIndex + 1; i < length; ++i) {
            tp1 = tempList.get(i);
            tp2 = tempList.get(i - 1);
            long t1 = tp1.getTimestamp();
            long t2 = tp2.getTimestamp();
            // if (distance(prePoint, tp1) < ((t1-2*kpTime+preTime)*(t1-2*kpTime+preTime) * SMAX * SMAX)) {
            //   break;
            // }
            // xi and xi-1 satisfy speed constraint
            if (distance(tp1, tp2) <= ((t1 - t2) * (t1 - t2) * SMAX * SMAX)) {
                if (top[i - 1] == -1) {
                    top[i] = i - 1;
                    len[i - 1]++;
                } else if (top[i - 1] > 0) {
                    top[i] = top[i - 1];
                    len[top[i - 1]]++;
                }
            } else { // xi vialote the speed constraint with xi-1
                for (int j = i - 1; j >= topIndex; j--) {
                    tp2 = tempList.get(j);
                    t2 = tp2.getTimestamp();
                    if ((distance(tp1, tp2) > ((t1 - t2) * (t1 - t2) * SMAX * SMAX) && top[j] > 0) || j == topIndex) {
                        if (distance(prePoint, tp1) <= ((t1 - preTime) * (t1 - preTime) * SMAX * SMAX)) {
                            // if (distance(prePoint, tp1) > ((t1-kpTime-1)*(t1-kpTime-1) * SMAX * SMAX)) {
                            // if (distance(prePoint, tp1) > ((t1-2*kpTime+preTime)*(t1-2*kpTime+preTime) * SMAX * SMAX)) {
                            //   // distance(prePoint, tp1) > ((t1-2*kpTime+preTime)*(t1-2*kpTime+preTime) * SMAX * SMAX)
                            //     top[i] = -1;
                            //     len[i] = 1;
                            // }
                            // else {
                            //     flag = true;
                            // }
                            top[i] = -1;
                            len[i] = 1;
                        }
                        break;
                    } else if (distance(tp1, tp2) > ((t1 - t2) * (t1 - t2) * SMAX * SMAX) && (top[j] == 0 || top[j] == -1)) {
                        continue;
                    } else if (distance(tp1, tp2) <= ((t1 - t2) * (t1 - t2) * SMAX * SMAX)) {
                        if (top[j] == -1) {
                            top[i] = j;
                            len[j]++;
                        } else if (top[j] > 0) {
                            top[i] = top[j];
                            len[top[j]]++;
                        }
                        break;
                    } else {
                        System.out.println("My2.local() error");
                        return;
                    }
                }
            }
            if (flag) {
                break;
            }
        }

        // find maxpoint
        int maxIndex = topIndex;
        for (int i = maxIndex; i < length; ++i) {
            if (len[i] > len[maxIndex]) {
                maxIndex = i;
            }
        }
        TimePoint maxPoint = tempList.get(maxIndex);
        long maxTime = maxPoint.getTimestamp();
//        double[] maxVal = new double[2];
//        maxVal = maxPoint.getModify();
        DMatrixRMaj maxVal = maxPoint.getValRepaired();

        // whether to repair
        boolean judge = judgeModify(preVal, maxVal, kpVal, preTime, maxTime, kpTime);

        // repair
        if (judge) {
//            double[] xy_modify = new double[2];
            DMatrixRMaj xy_modify = new DMatrixRMaj(dimension, 1);
            double pre_dis = this.SMAX * (kpTime - preTime);
            // double next_dis = this.SMAX * (maxTime-kpTime);
            // double pre_next_dis = calculateDistance(maxPoint, prePoint);
            double pre_next_dis = this.SMAX * (maxTime - preTime);
            // double pre_kp_dis = calculateDistance(kp, prePoint);
            if (pre_next_dis > pre_dis) {
                double tmp1 = (maxVal.get(0,0) - preVal.get(0,0)) * (pre_dis / pre_next_dis) + preVal.get(0,0);
                double tmp2 = (maxVal.get(1,0) - preVal.get(1,0)) * (pre_dis / pre_next_dis) + preVal.get(1,0);
                xy_modify.set(0,0, tmp1);
                xy_modify.set(1,0, tmp2);
            } else {
                // else if(pre_next_dis > pre_dis && pre_next_dis > next_dis){
                double tmp1 = (maxVal.get(0,0) + preVal.get(0,0)) / 2;
                double tmp2 = (maxVal.get(1,0) + preVal.get(1,0)) / 2;
                xy_modify.set(0,0, tmp1);
                xy_modify.set(1,0, tmp2);

            }
//            kp.setModify(xy_modify);
            kp.setValRepaired(new DMatrixRMaj(xy_modify));
        }
    }
}
