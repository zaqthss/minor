package com.MINOR.entity;

@Deprecated
public class Discord {
    public int index;
    public int length;
    public double d;

    public Discord(int index, int length, double d) {
        this.index = index;
        this.length = length;
        this.d = d;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[" + index + "," + (index + length) + "]" + ",L=" + length + ",d=" + d);
        return sb.toString();
    }


}
