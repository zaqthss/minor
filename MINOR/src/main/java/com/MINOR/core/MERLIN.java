package com.MINOR.core;

import com.MINOR.Utils.MyMathUtils;
import com.MINOR.entity.Discord;
import com.MINOR.entity.TimeSeries;

import java.util.*;
import java.util.concurrent.CopyOnWriteArraySet;

@Deprecated
public class MERLIN {
    private Set<Integer> C;
    private Set<Discord> D;
    private TimeSeries ts;      // 正则化后的时序数据
    private int tsLength;
    private int dim;

    public MERLIN(TimeSeries ts) {
        this.ts = ts;

        tsLength = ts.size();
        dim = ts.getDimension();
    }

    private void run(int minL, int maxL) {
        double r = 2 * Math.sqrt(minL);
        double d_minL = Double.NEGATIVE_INFINITY;
        while (d_minL < 0) {
            DRAG(minL, r);
            r /= 2;
        }
        for (int i = minL + 1; i < minL + 4; i++) {
            double di = Double.NEGATIVE_INFINITY;
            while (di < 0) {
                r = 0.99 * di;
                r *= 0.99;
            }
        }

    }

    public void DRAG(int L, double r) {
        getCandidates(L, r);
        discordsRefinement(L, r);
        for (Discord discord : D) {
            for (int i = discord.index; i < discord.index + discord.length; i++) {
                ts.setAnomalyAt(i, true);
            }
        }
    }

    private boolean getCandidates(int L, double r) {
        C = new CopyOnWriteArraySet<>();
        boolean isCandidate;
        double d;
        TimeSeries tsi, tsj;
        for (int i = 0; i < tsLength - L; i ++) {
            isCandidate = true;
            if (!C.isEmpty()) {
                for (int j : C) {
                    if (Math.abs(i - j) >= L) {
                        tsi = ts.getNormSubList(i, L);
                        // non-self match
                        tsj = ts.getNormSubList(j, L);
                        d = MyMathUtils.eulerDistObs(tsi, tsj);
                        if (d < r) {
                            C.remove(j);
                            isCandidate = false;
                        }
                    }else{
                        isCandidate = false;
                    }
                }
            }
            if (isCandidate) {
                C.add(i);
            }
        }
        if (C.isEmpty()) {
            return false;
        } else {
            return true;
        }
    }

    private void discordsRefinement(int L, double r) {
        D = new CopyOnWriteArraySet<>();
        boolean isDiscord;
        Map<Integer, Boolean> isD = new HashMap<>();
        double d, dj;
        TimeSeries tsi, tsj;
        if (!C.isEmpty()) {
            for (int j : C) {
                dj = Double.MAX_VALUE;
                isDiscord = true;
                tsj = ts.getNormSubList(j, L);
                for (int i = 0; i < tsLength - L; i++) {
                    if (Math.abs(i - j) >= L) {
                        tsi = ts.getNormSubList(i, L);
                        d = MyMathUtils.eulerDistObs(tsi, tsj);
                        if (d < r) {
                            C.remove(j);
                            isDiscord = false;
                        } else {
                            dj = Math.min(dj, d);
                        }
                    }
                }
                if (isDiscord) {
                    D.add(new Discord(j, L, dj));
                }
            }
        }
//        for (int i = 0; i < tsLength - L; i += L) {
//            isDiscord = true;
//            di = Double.MAX_VALUE;
//            if (!C.isEmpty()) {
//                for (int j : C) {
//                    if (Math.abs(i - j) >= L) {
//                        tsi = ts.getNormSubList(i, L);
//                        tsj = ts.getNormSubList(j, L);
//                        d = CustomMathUtils.eulerDistObs(tsi, tsj);
//                        if (d < r) {
//                            C.remove(j);
//                            isDiscord = false;
//                        } else {
//                            di = Math.min(di, d);
//                        }
//                    }
//                }
//            }
//            if (isDiscord) {
//                D.add(new Discord(i, L, di));
//            }
//        }
    }

    /**
     * 将D中的Discord按照index升序打印
     */
    public void printDiscordSet() {
        List<Discord> lD = new ArrayList<>(D);
        Collections.sort(lD, Comparator.comparingInt(d -> d.index));
        if(lD.isEmpty()){
            System.out.println("[INFO]: No discord was found.");
            return;
        }
        for (Discord d : lD) {
            System.out.println(d);
        }
    }
}
