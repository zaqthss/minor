package com.MINOR.entity;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * 经过排序的数组，可获取指定百分位的数据
 */
public class SortedList<T> {
    private List<T> list;
    private Comparator<? super T> comparator;

    public SortedList(Comparator<? super T> comparator) {
        this.list = new ArrayList<T>();
        this.comparator = comparator;
    }

    public void add(T value) {
        int index = Collections.binarySearch(list, value, comparator);
        if (index < 0) {
            // 如果元素不存在，binarySearch 会返回 (-插入点 - 1)
            index = -index - 1;
        }
        list.add(index, value);
    }

    public List<T> getList() {
        return list;
    }

    public T getPercentileVal(int percentile) {
        int index = (int) Math.ceil((double) percentile / 100 * list.size()) - 1;
        index = Math.max(0, index);
        index = Math.min(list.size() - 1, index);
        return list.get(index);
    }
}
