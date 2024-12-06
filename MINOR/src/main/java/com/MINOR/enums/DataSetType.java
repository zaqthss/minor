package com.MINOR.enums;

public enum DataSetType {
    ILD("ILD"),
    ECG("ECG");

    private String value;

    private DataSetType(String value) {
        this.value = value;
    }

    @Override
    public String toString() {
        return value;
    }
}
