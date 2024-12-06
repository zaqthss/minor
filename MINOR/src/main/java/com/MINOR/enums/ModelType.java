package com.MINOR.enums;

public enum ModelType {
    MD_AR("AR"),
    MD_ARX("ARX"),
    MD_IMR("IMR"),
    MD_VAR("VAR"),
    MD_VARX("VARX"),
    MD_MINOR_U("MINOR-U"),
    MD_MINOR_B("MINOR-B"),
    MD_MINOR_BUNI("MINOR-B-Uni"),
    MD_MINOR_O("MINOR-O");

    private final String displayName;

    ModelType(String displayName) {
        this.displayName = displayName;
    }

    @Override
    public String toString() {
        return this.displayName;
    }
}
