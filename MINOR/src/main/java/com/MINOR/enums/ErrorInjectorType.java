package com.MINOR.enums;

/**
 * The error injection type specified during the construction of the ErrorInjection class.
 */
public enum ErrorInjectorType {
    SHIFT("SHIFT"),
    INNOVATIONAL("INNO");

    private final String value;

    ErrorInjectorType(String value) {
        this.value = value;
    }

    @Override
    public String toString() {
        return value;
    }
}
