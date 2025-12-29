package org.igv.sam;

import htsjdk.tribble.Feature;

/**
 * Genomic interval showing which areas have had reads removed (downsampled)
* @author jrobinso
*         Date: 8/16/12
*         Time: 10:03 PM
*/
public class DownsampledInterval implements Feature {
    private int start;
    private int end;
    private int count;

    public DownsampledInterval(int start, int end, int count) {
        this.start = start;
        this.end = end;
        this.count = count;
    }

    public String toString() {
        return start + "-" + end + " (" + count + ")";
    }

    public int getCount() {
        return count;
    }

    public int getEnd() {
        return end;
    }

    public int getStart() {
        return start;
    }

    public String getChr() {
        return null;
    }

    @Override
    public String getContig() {
        return null;
    }

    public String getValueString() {
        return "Interval [" + start + "-" + end + "] <br>" + count + " reads removed.";
    }

    public void incCount() {
        this.count += 1;
    }
}
