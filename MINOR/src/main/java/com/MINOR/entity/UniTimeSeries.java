package com.MINOR.entity;

import com.MINOR.Utils.MyMathUtils;

import java.util.ArrayList;
import java.util.List;

public class UniTimeSeries {
    public List<UniTimePoint> s;
    public List<Integer> candidateList;

    public UniTimeSeries(List<UniTimePoint> s) {
        this.s = s;
        candidateList = new ArrayList<>();
    }

    public UniTimeSeries() {
        this.s = new ArrayList<>();
        candidateList = new ArrayList<>();
    }

    public void add(UniTimePoint p) {
        s.add(p);
    }

    // assert that the length matched
    public static UniTimeSeries getFromList(List<Long> timestamp, List<Double> obs, List<Double> truth, List<Integer> labels) {
        int size = timestamp.size();
        List<UniTimePoint> ss = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            UniTimePoint tp = new UniTimePoint(timestamp.get(i), truth.get(i), obs.get(i), Double.MAX_VALUE, obs.get(i), MyMathUtils.int2Bool(labels.get(i)));
            ss.add(tp);
        }
        return new UniTimeSeries(ss);
    }

    public List<Double> toList() {
        List<Double> result = new ArrayList<>();
        for (UniTimePoint uniTimePoint : s) {
            result.add(uniTimePoint.getValRepaired());
        }
        return result;
    }

    /**
     * 将candidateList中所有数据点的修复候选值清除，并清空candidateList
     */
    public void clearCandidates() {
        for (Integer i : candidateList) {
            s.get(i).setValCandidate(Double.MAX_VALUE);
        }
        candidateList.clear();
    }

    public int size() {
        return s.size();
    }

    public void concat(UniTimeSeries _ts) {
        this.s.addAll(_ts.s);
    }

    /**
     * 获得前n个数据组成的子序列，自身会去除该部分
     *
     * @param n 子序列长度
     */
    public UniTimeSeries popSubList(int n) {
        List<UniTimePoint> subS = s.subList(0, n);
        this.s = s.subList(n, s.size());
        return new UniTimeSeries(subS);
    }

    /**
     * 将原序列分割成若干子序列。每个子序列可以以零个、一个或多个连续的未标记数据开头,随后紧跟着若干已标记数据。
     */
    public List<UniTimeSeries> breakdownByLabel() {
        List<UniTimeSeries> result = new ArrayList<>();
        List<UniTimePoint> currentSegment = new ArrayList<>();
        boolean foundLabeledPoint = false;
        boolean notLabeled;
        for (UniTimePoint point : s) {
            notLabeled = !point.getLabel();
            if (notLabeled) {
                // for unlabeled points
                if (foundLabeledPoint) {
                    // form a new sub list
                    foundLabeledPoint = false;
                    result.add(new UniTimeSeries(currentSegment));
                    currentSegment = new ArrayList<>();
                }
            } else {
                // for labeled points
                foundLabeledPoint = true;
            }
            currentSegment.add(point);
        }
        if (!currentSegment.isEmpty()) {
            result.add(new UniTimeSeries(currentSegment));
        }
        return result;
    }

    /**
     * 分解为长度为l的子序列
     */
    public List<UniTimeSeries> breakdown(int l) {
        List<UniTimePoint> subS;
        List<UniTimeSeries> ls = new ArrayList<>();
        int n = s.size();
        int i;
        for (i = 0; i < n; i += l) {
            subS = new ArrayList<>();
            for (int j = 0; j < l; j++) {
                subS.add(s.get(i + j));
            }
            ls.add(new UniTimeSeries(subS));
        }
        subS = new ArrayList<>();
        for (i -= l; i < n; i++) {
            subS.add(s.get(i));
        }
        if (!subS.isEmpty()) {
            ls.add(new UniTimeSeries(subS));
        }
        return ls;
    }

    public UniTimeSeries reverse() {
        int n = s.size();
        List<UniTimePoint> rs = new ArrayList<>();
        for (int i = n - 1; i >= 0; i--) {
            rs.add(new UniTimePoint(s.get(i)));
        }
        return new UniTimeSeries(rs);
    }

    public UniTimePoint getP(int index) {
        return s.get(index);
    }

    public long getTimestampAt(int index) {
        return s.get(index).getTimestamp();
    }

    public double getTruthAt(int index) {
        return s.get(index).getValTruth();
    }

    public void setTruthAt(int index, double val) {
        s.get(index).setValTruth(val);
    }

    public double getObsAt(int index) {
        return s.get(index).getValObs();
    }

    public void setObsAt(int index, double val) {
        s.get(index).setValObs(val);
    }

    public void setCandidateAt(int index, double val) {
        candidateList.add(index);
        s.get(index).setValCandidate(val);
    }

    public double getCandidateAt(int index) {
        return s.get(index).getValCandidate();
    }

    public void setRepairAt(int index, double val) {
        UniTimePoint p = s.get(index);
        if (Double.MAX_VALUE == p.getValRepaired()) {
            p.setValLastModify(val);
        } else {
            p.setValLastModify(p.getValRepaired());
        }
        p.setValRepaired(val);
    }

    public double getRepairAt(int index) {
        return s.get(index).getValRepaired();
    }

    public void setLastRepairAt(int index, double val) {
        s.get(index).setValLastRepaired(val);
    }

    public double getLastRepairAt(int index) {
        return s.get(index).getValLastRepaired();
    }

    public Boolean getLabelAt(int index) {
        return s.get(index).getLabel();
    }

    public void setLabelAt(int index, boolean val) {
        s.get(index).setLabel(val);
    }

    public double getResidualAt(int index) {
        return s.get(index).getValRepaired() - s.get(index).getValObs();
    }

    public double getLastResidualAt(int index) {
        return s.get(index).getValLastRepaired() - s.get(index).getValObs();
    }

    public double getLastModifyAt(int index) {
        return s.get(index).getValLastModify();
    }

    public void setLastModifyAt(int index, double val) {
        s.get(index).setValLastModify(val);
    }
}
