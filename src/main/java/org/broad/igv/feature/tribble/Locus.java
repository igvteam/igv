package org.broad.igv.feature.tribble;

import htsjdk.tribble.Feature;

/**
 * A minimal representation of a Feature.
 */
public class Locus implements Feature {

    String chr;
    int start;
    int end;

    /**
     *
     * @param chr
     * @param start
     * @param end
     */
    public Locus(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

    public String getChr() {
        return chr;
    }

    @Override
    public String getContig() {
        return chr;
    }

    public int getStart() {
        return start;  
    }

    public int getEnd() {
        return end;
    }
}
