package org.broad.igv.ucsc.bb;

import org.broad.igv.feature.LocusScore;

public class WigDatum implements LocusScore {

    private String chr;
    private int start;
    private int end;
    private float score;

    public WigDatum(String chr, int start, int end, float score) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.score = score;
    }

    @Override
    public float getScore() {
        return score;
    }

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    @Override
    public String getChr() {
        return chr;
    }
}
