package com.MINOR.entity;

import com.MINOR.Utils.MyMathUtils;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

import java.util.ArrayList;
import java.util.List;

public class TimeSeries {
    public List<TimePoint> s;
    public List<Integer> candidateList;

    public TimeSeries(List<TimePoint> s) {
        this.s = s;
        candidateList = new ArrayList<>();
    }

    public TimeSeries() {
        this.s = new ArrayList<>();
        candidateList = new ArrayList<>();
    }

    public void add(TimePoint p) {
        s.add(p);
    }

    public void clearRepair(){
        for(TimePoint tp : s){
            tp.setValCandidate(null);
            tp.setValLastRepaired(new DMatrixRMaj(tp.getValObs()));
            tp.setValLastModify(new DMatrixRMaj(tp.getValObs()));
            if(!tp.getLabel()){
                tp.setValRepaired(new DMatrixRMaj(tp.getValObs()));
            }
        }
    }

    // assert that the length matched
    public static TimeSeries getFromList(List<Long> timestamp, List<Double> obs, List<Double> truth, List<Integer> labels) {
        int size = timestamp.size();
        List<TimePoint> ss = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            double[] obsMat = {obs.get(i)};
            double[] truthMat = {truth.get(i)};
            TimePoint tp = new TimePoint(timestamp.get(i), 1, new DMatrixRMaj(truthMat), new DMatrixRMaj(obsMat), null, null, MyMathUtils.int2Bool(labels.get(i)));
            ss.add(tp);
        }
        return new TimeSeries(ss);
    }

    public List<Double> toList() {
        if (this.getDimension() != 1) {
            System.out.println("[ERROR] dimension must be 1.");
            return null;
        }
        List<Double> result = new ArrayList<>();
        for (TimePoint timePoint : s) {
            result.add(timePoint.getValRepaired().get(0, 0));
        }
        return result;
    }

    /**
     * turn a multi-variate time series to n uni-variate time series
     */
    public List<UniTimeSeries> toUniTimeSeries() {
        int dim = getDimension();
        int size = s.size();
        List<UniTimeSeries> unis = new ArrayList<>();
        for (int j = 0; j < dim; j++) {
            unis.add(new UniTimeSeries());
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < dim; j++) {
                double c = this.getCandidateAt(i) == null ? Double.MAX_VALUE : this.getCandidateAt(i).get(j, 0);
                double r = this.getRepairAt(i) == null ? Double.MAX_VALUE : this.getRepairAt(i).get(j, 0);
                UniTimePoint p = new UniTimePoint(this.getTimestampAt(i),
                        this.getTruthAt(i).get(j, 0),
                        this.getObsAt(i).get(j, 0),
                        c,
                        r,
                        this.getLabelAt(i));
                unis.get(j).s.add(p);
            }
        }
        return unis;
    }

    /**
     * turn a multi-variate time series to n uni-variate time series
     */
    public List<TimeSeries> toSinglePointTimeSeries() {
        int dim = getDimension();
        int size = s.size();
        List<TimeSeries> unis = new ArrayList<>();
        for (int j = 0; j < dim; j++) {
            unis.add(new TimeSeries());
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < dim; j++) {
                DMatrixRMaj c = null;
                DMatrixRMaj r = null;
                double[] t_ = new double[1];
                t_[0] = this.getTruthAt(i).get(j, 0);
                double[] j_ = new double[1];
                j_[0] = this.getObsAt(i).get(j, 0);
                TimePoint p = new TimePoint(this.getTimestampAt(i),
                        1,
                        new DMatrixRMaj(t_),
                        new DMatrixRMaj(j_),
                        c,
                        r,
                        this.getLabelAt(i));
                unis.get(j).s.add(p);
            }
        }
        return unis;
    }

    /**
     * turn some uni-variate time series to a multi-variate time series
     */
    public static TimeSeries aggregate(List<UniTimeSeries> unis) {
        int size = unis.get(0).s.size();
        int dim = unis.size();
        TimeSeries ts = new TimeSeries();
        boolean b = false;
        for (int i = 0; i < size; i++) {
            long timestamps = unis.get(0).s.get(i).getTimestamp();
            DMatrixRMaj vt = new DMatrixRMaj(dim, 1);
            DMatrixRMaj vo = new DMatrixRMaj(dim, 1);
            DMatrixRMaj vc = new DMatrixRMaj(dim, 1);
            DMatrixRMaj vr = new DMatrixRMaj(dim, 1);
            for (int j = 0; j < dim; j++) {
                UniTimePoint up = unis.get(j).s.get(i);
                vt.set(j, 0, up.getValTruth());
                vo.set(j, 0, up.getValObs());
                vc.set(j, 0, up.getValCandidate() == Double.MAX_VALUE ? up.getValObs() : up.getValCandidate());
                vr.set(j, 0, up.getValRepaired() == Double.MAX_VALUE ? up.getValObs() : up.getValRepaired());
                b = up.getLabel();
            }
            TimePoint p = new TimePoint(timestamps, dim, vt, vo, vc, vr, b);
            ts.s.add(p);
        }
        return ts;
    }

    public static TimeSeries aggregateSingleTs(List<TimeSeries> unis) {
        int size = unis.get(0).s.size();
        int dim = unis.size();
        TimeSeries ts = new TimeSeries();
        boolean b = false;
        for (int i = 0; i < size; i++) {
            long timestamps = unis.get(0).s.get(i).getTimestamp();
            DMatrixRMaj vt = new DMatrixRMaj(dim, 1);
            DMatrixRMaj vo = new DMatrixRMaj(dim, 1);
            DMatrixRMaj vc = new DMatrixRMaj(dim, 1);
            DMatrixRMaj vr = new DMatrixRMaj(dim, 1);
            for (int j = 0; j < dim; j++) {
                TimePoint up = unis.get(j).s.get(i);
                vt.set(j, 0, up.getValTruth().get(0, 0));
                vo.set(j, 0, up.getValObs().get(0, 0));
                vc.set(j, 0, up.getValCandidate() == null ? up.getValObs().get(0, 0) : up.getValCandidate().get(0, 0));
                vr.set(j, 0, up.getValRepaired() == null ? up.getValObs().get(0, 0) : up.getValRepaired().get(0, 0));
                b = up.getLabel();
            }
            TimePoint p = new TimePoint(timestamps, dim, vt, vo, vc, vr, b);
            ts.s.add(p);
        }
        return ts;
    }

    /**
     * 检查序列中点的维度是否一致
     *
     * @return true if consistent
     */
    public boolean checkDimension() {
        return s.stream()
                .map(TimePoint::getDimension)
                .distinct()
                .count() == 1;
    }

    /**
     * 将candidateList中所有数据点的修复候选值清除，并清空candidateList
     */
    public void clearCandidates() {
        for (Integer i : candidateList) {
            s.get(i).setValCandidate(null);
        }
        candidateList.clear();
    }

    public int size() {
        return s.size();
    }

    public void concat(TimeSeries _ts) {
        this.s.addAll(_ts.s);
    }

    /**
     * 获得前n个数据组成的子序列，自身会去除该部分
     *
     * @param n 子序列长度
     * @return 子序列实例对象
     */
    public TimeSeries popSubList(int n) {
        List<TimePoint> subS = s.subList(0, n);
        this.s = s.subList(n, s.size());
        return new TimeSeries(subS);
    }

    /**
     * 获得长度为L的子序列实例
     *
     * @param start 起点序号
     * @param L     子序列长度
     * @return 子序列实例对象
     */
    public TimeSeries getSubList(int start, int L) {
        List<TimePoint> list = s.subList(start, start + L);
        return new TimeSeries(list);
    }

    /**
     * 获得长度为L并且标准化的子序列实例
     *
     * @param start 起点序号
     * @param L     子序列长度
     * @return 标准化的子序列实例对象
     */
    public TimeSeries getNormSubList(int start, int L) {
        List<TimePoint> list = s.subList(start, start + L);
        TimeSeries ts2 = new TimeSeries(list);
        return ts2.getNormalizedTS();
    }

    /**
     * 将原序列分割成若干子序列。每个子序列可以以零个、一个或多个连续的未标记数据开头,随后紧跟着若干已标记数据。
     */
    public List<TimeSeries> breakdownByLabel() {
        List<TimeSeries> result = new ArrayList<>();
        List<TimePoint> currentSegment = new ArrayList<>();
        boolean foundLabeledPoint = false;
        for (TimePoint point : s) {
            if (!point.getLabel()) {
                // for unlabeled points
                if (foundLabeledPoint) {
                    // form a new sub list
                    foundLabeledPoint = false;
                    result.add(new TimeSeries(currentSegment));
                    currentSegment = new ArrayList<>();
                }
            } else {
                // for labeled points
                foundLabeledPoint = true;
            }
            currentSegment.add(point);
        }
        if (!currentSegment.isEmpty()) {
            result.add(new TimeSeries(currentSegment));
        }
        return result;
    }

    /**
     * 分解为长度为1的子序列
     */
    public List<TimeSeries> breakdown() {
        List<TimePoint> subS;
        List<TimeSeries> ls = new ArrayList<>();
        for (TimePoint timePoint : s) {
            subS = new ArrayList<>();
            subS.add(timePoint);
            ls.add(new TimeSeries(subS));
        }
        return ls;
    }

    /**
     * 分解为长度为l的子序列
     */
    public List<TimeSeries> breakdown(int l) {
        List<TimePoint> subS;
        List<TimeSeries> ls = new ArrayList<>();
        int n = s.size();
        int i;
        for (i = 0; i < n; i += l) {
            subS = new ArrayList<>();
            for (int j = 0; j < l; j++) {
                subS.add(s.get(i + j));
            }
            ls.add(new TimeSeries(subS));
        }
        subS = new ArrayList<>();
        for (i -= l; i < n; i++) {
            subS.add(s.get(i));
        }
        if (!subS.isEmpty()) {
            ls.add(new TimeSeries(subS));
        }
        return ls;
    }

    /**
     * 获得一个TimePoint list取反后的新TimeSeries对象，
     */
    public TimeSeries reverse() {
        int n = s.size();
        List<TimePoint> rs = new ArrayList<>();
        for (int i = n - 1; i >= 0; i--) {
            rs.add(new TimePoint(s.get(i)));
        }
        return new TimeSeries(rs);
    }

    /**
     * 对时序数据进行标准化
     *
     * @return 标准化后的新时序数据实例
     */
    public TimeSeries getNormalizedTS() {
        List<TimePoint> ls = new ArrayList<>();
        TimePoint tp;
        int n = size();
        int m = getDimension();
        double[] sum = new double[m];
        double[] miu = new double[m];
        double[] sigma = new double[m];
        DMatrixRMaj x = new DMatrixRMaj(m, 1);
        double val;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                val = getObsAt(i).get(j, 0);
                sum[j] += val;
            }
        }
        for (int i = 0; i < m; i++) {
            // 期望
            miu[i] = sum[i] / n;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                val = getObsAt(i).get(j, 0);
                sigma[j] += Math.pow(val - miu[j], 2);
            }
        }
        for (int i = 0; i < m; i++) {
            // 标准差
            sigma[i] = Math.sqrt(sigma[i] / n);
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                val = (getObsAt(i).get(j, 0) - miu[j]) / sigma[j];
                x.set(j, 0, val);
            }
            tp = new TimePoint(s.get(i));
            tp.setValObs(new DMatrixRMaj(x));
            ls.add(i, tp);
        }
        return new TimeSeries(ls);
    }

    /**
     * 应用anomalyDetection的结果，将detectedAsNormal赋值给label，如果原本label就为true，则置detectedAsNormal为false
     */
    public void applyAD() {
        for (TimePoint p : this.s) {
            if (p.isDetectedAsNormal()) {
                if (!p.getLabel()) {
                    p.setLabel(true);
                } else {
                    p.setDetectedAsNormal(false);
                }
            }
        }
    }

    /**
     * 撤销anomalyDetection的结果，将detectedAsNormal为true者，label置为false
     */
    public void revertAD() {
        for (TimePoint p : this.s) {
            if (p.isDetectedAsNormal()) {
                p.setDetectedAsNormal(false);
                p.setLabel(false);
            }
        }
    }

    public int getDimension() {
        return s.get(0).getDimension();
    }

    public TimePoint getP(int index) {
        return s.get(index);
    }

    public long getTimestampAt(int index) {
        return s.get(index).getTimestamp();
    }

    public DMatrixRMaj getTruthAt(int index) {
        return s.get(index).getValTruth();
    }

    public DMatrixRMaj getObsAt(int index) {
        return s.get(index).getValObs();
    }

    public void setCandidateAt(int index, DMatrixRMaj mat) {
        candidateList.add(index);
        s.get(index).setValCandidate(mat);
    }

    public DMatrixRMaj getCandidateAt(int index) {
        return s.get(index).getValCandidate();
    }

    /**
     * 设置修复值，会自动记录修复之前的值并存储到valLastModify
     *
     * @param index 下标，从0开始
     * @param mat   数据点的值
     */
    public void setRepairAt(int index, DMatrixRMaj mat) {
        TimePoint p = s.get(index);
        DMatrixRMaj vec = p.getValRepaired();
        if (null == vec) {
            p.setValLastModify(new DMatrixRMaj(mat));
        } else {
            p.setValLastModify(new DMatrixRMaj(vec));
        }
        p.setValRepaired(mat);
    }

    public DMatrixRMaj getRepairAt(int index) {
        return s.get(index).getValRepaired();
    }

    public void setLastRepairAt(int index, DMatrixRMaj mat) {
        s.get(index).setValLastRepaired(mat);
    }

    public DMatrixRMaj getLastRepairAt(int index) {
        return s.get(index).getValLastRepaired();
    }

    public void setLastModifyAt(int index, DMatrixRMaj mat) {
        s.get(index).setValLastModify(mat);
    }

    public DMatrixRMaj getLastModifyAt(int index) {
        return s.get(index).getValLastModify();
    }

    public void setLabelAt(int index, boolean lb) {
        s.get(index).setLabel(lb);
    }

    public boolean getLabelAt(int index) {
        return s.get(index).getLabel();
    }

    public DMatrixRMaj getResidualAt(int index) {
        DMatrixRMaj Z = new DMatrixRMaj(s.get(index).getValObs().numRows, s.get(index).getValObs().numCols);
        return CommonOps_DDRM.subtract(s.get(index).getValRepaired(), s.get(index).getValObs(), Z);
    }

    public DMatrixRMaj getLastResidualAt(int index) {
        DMatrixRMaj Z = new DMatrixRMaj(s.get(index).getValObs().numRows, s.get(index).getValObs().numCols);
        return CommonOps_DDRM.subtract(s.get(index).getValLastRepaired(), s.get(index).getValObs(), Z);
    }

    public boolean getAnomalyAt(int index) {
        return s.get(index).anomaly;
    }

    public void setAnomalyAt(int index, boolean b) {
        s.get(index).anomaly = b;
    }

    public boolean isDetectedAsNormal(int index) {
        return s.get(index).isDetectedAsNormal();
    }

    public void setDetectedAsNormal(int index, boolean b) {
        s.get(index).setDetectedAsNormal(b);
    }
}
