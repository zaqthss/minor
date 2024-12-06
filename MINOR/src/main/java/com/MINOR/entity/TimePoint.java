package com.MINOR.entity;

import org.ejml.data.DMatrixRMaj;

public class TimePoint {
    private long timestamp;
    private int dimension;
    private DMatrixRMaj valTruth;
    private DMatrixRMaj valObs;
    private DMatrixRMaj valCandidate;
    private DMatrixRMaj valRepaired;
    /**
     * 这个值用于MIMR基于某一个点上一轮迭代的值，在增量式参数估计中使用
     */
    private DMatrixRMaj valLastRepaired;
    /**
     * 这个值用于MIMR_bipv记录某一个点在上一次修改后的值，用于计算两次修改向量的夹角
     */
    private DMatrixRMaj valLastModify;
    private boolean label;
    private boolean detectedAsNormal;
    public boolean anomaly;

    /**
     * 一个m维的时序数据点，应当以列向量形式表示
     *
     * @param timestamp    timestamp
     * @param m            dimension
     * @param valTruth     ground truth
     * @param valObs       observed value
     * @param valCandidate candidate for repairing
     * @param valRepaired  repaired value
     * @param label        true if labeled
     */
    public TimePoint(long timestamp, int m, DMatrixRMaj valTruth, DMatrixRMaj valObs, DMatrixRMaj valCandidate, DMatrixRMaj valRepaired, boolean label) {
        this.timestamp = timestamp;
        this.dimension = m;
        this.valTruth = valTruth;
        this.valObs = valObs;
        this.valCandidate = valCandidate;
        this.valRepaired = valRepaired;
        this.valLastRepaired = valRepaired;
        this.valLastModify = valRepaired;
        this.label = label;
        this.detectedAsNormal = false;
        anomaly = false;
    }

    public TimePoint(TimePoint p) {
        this.timestamp = p.timestamp;
        this.dimension = p.dimension;
        this.valTruth = new DMatrixRMaj(p.valTruth);
        this.valObs = new DMatrixRMaj(p.valObs);
        if (p.valCandidate != null) {
            this.valCandidate = new DMatrixRMaj(p.valCandidate);
        } else {
            this.valCandidate = null;
        }
        if (p.valRepaired != null) {
            this.valRepaired = new DMatrixRMaj(p.valRepaired);
        } else {
            this.valRepaired = null;
        }
        if (p.valLastRepaired != null) {
            this.valLastRepaired = new DMatrixRMaj(p.valLastRepaired);
        } else {
            if(p.valRepaired != null){
                this.valLastRepaired = valRepaired;
            }else{
                this.valLastRepaired = null;
            }
        }
        if (p.valLastModify != null) {
            this.valLastModify = new DMatrixRMaj(p.valLastModify);
        } else {
            if(p.valRepaired != null){
                this.valLastModify = valRepaired;
            }else{
                this.valLastModify = null;
            }
        }
        this.label = p.label;
        this.detectedAsNormal = p.detectedAsNormal;
        anomaly = false;
    }

    public long getTimestamp() {
        return timestamp;
    }

    public void setTimestamp(long timestamp) {
        this.timestamp = timestamp;
    }

    public int getDimension() {
        return dimension;
    }

    public void setDimension(int dimension) {
        this.dimension = dimension;
    }

    public DMatrixRMaj getValTruth() {
        return valTruth;
    }

    public void setValTruth(DMatrixRMaj valTruth) {
        this.valTruth = valTruth;
    }

    public DMatrixRMaj getValObs() {
        return valObs;
    }

    public void setValObs(DMatrixRMaj valObs) {
        this.valObs = valObs;
    }

    public DMatrixRMaj getValCandidate() {
        return valCandidate;
    }

    public void setValCandidate(DMatrixRMaj valCandidate) {
        this.valCandidate = valCandidate;
    }

    public DMatrixRMaj getValRepaired() {
        return valRepaired;
    }

    public void setValRepaired(DMatrixRMaj valRepaired) {
        this.valRepaired = valRepaired;
    }

    public boolean getLabel() {
        return label;
    }

    public void setLabel(boolean label) {
        this.label = label;
    }

    public DMatrixRMaj getValLastRepaired() {
        return valLastRepaired;
    }

    public void setValLastRepaired(DMatrixRMaj valLastRepaired) {
        this.valLastRepaired = valLastRepaired;
    }

    public DMatrixRMaj getValLastModify() {
        return valLastModify;
    }

    public void setValLastModify(DMatrixRMaj valLastModify) {
        this.valLastModify = valLastModify;
    }

    public boolean isDetectedAsNormal() {
        return detectedAsNormal;
    }

    public void setDetectedAsNormal(boolean detectedAsNormal) {
        this.detectedAsNormal = detectedAsNormal;
    }
}
