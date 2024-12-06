package com.MINOR.entity;

public class UniTimePoint {
    private long timestamp;
    private double valTruth;
    private double valObs;
    private double valCandidate;
    private double valRepaired;
    private double valLastRepaired;
    private double valLastModify;
    private Boolean label;

    /**
     * 用Double.MAX_VALUE表示不存在
     *
     * @param timestamp    timestamp
     * @param valTruth     ground truth
     * @param valObs       observed value
     * @param valCandidate candidate for repairing
     * @param valRepaired  repaired value
     * @param label        true if labeled
     */
    public UniTimePoint(long timestamp, double valTruth, double valObs, double valCandidate, double valRepaired, Boolean label) {
        this.timestamp = timestamp;
        this.valTruth = valTruth;
        this.valObs = valObs;
        this.valCandidate = valCandidate;
        this.valRepaired = valRepaired;
        this.valLastRepaired = valRepaired;
        this.valLastModify = valRepaired;
        this.label = label;
    }

    /**
     * 用Double.MAX_VALUE表示不存在
     *
     * @param timestamp    timestamp
     * @param valObs       observed value
     */
    public UniTimePoint(long timestamp, double valObs) {
        this.timestamp = timestamp;
        this.valTruth = Double.MAX_VALUE;
        this.valObs = valObs;
        this.valCandidate = Double.MAX_VALUE;
        this.valRepaired = Double.MAX_VALUE;
        this.valLastRepaired = valRepaired;
        this.valLastModify = valRepaired;
        this.label = false;
    }

    public UniTimePoint(UniTimePoint p) {
        this.timestamp = p.timestamp;
        this.valTruth = p.valTruth;
        this.valObs = p.valObs;
        this.valCandidate = p.valCandidate;
        this.valRepaired = p.valRepaired;
        this.valLastRepaired = p.valLastRepaired;
        this.valLastModify = p.valLastModify;
        this.label = p.label;
    }

    public long getTimestamp() {
        return timestamp;
    }

    public void setTimestamp(long timestamp) {
        this.timestamp = timestamp;
    }

    public double getValTruth() {
        return valTruth;
    }

    public void setValTruth(double valTruth) {
        this.valTruth = valTruth;
    }

    public double getValObs() {
        return valObs;
    }

    public void setValObs(double valObs) {
        this.valObs = valObs;
    }

    public double getValCandidate() {
        return valCandidate;
    }

    public void setValCandidate(double valCandidate) {
        this.valCandidate = valCandidate;
    }

    public double getValRepaired() {
        return valRepaired;
    }

    public void setValRepaired(double valRepaired) {
        this.valRepaired = valRepaired;
    }

    public double getValLastRepaired() {
        return valLastRepaired;
    }

    public void setValLastRepaired(double valLastRepaired) {
        this.valLastRepaired = valLastRepaired;
    }

    public Boolean getLabel() {
        return label;
    }

    public void setLabel(Boolean label) {
        this.label = label;
    }

    public double getValLastModify() {
        return valLastModify;
    }

    public void setValLastModify(double valLastModify) {
        this.valLastModify = valLastModify;
    }
}
