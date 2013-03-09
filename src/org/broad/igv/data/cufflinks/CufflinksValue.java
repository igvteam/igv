package org.broad.igv.data.cufflinks;

import org.broad.igv.feature.LocusScore;

/**
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 9:32 PM
 */
abstract public class CufflinksValue implements LocusScore {
    String chr;
    int start;
    int end;

    public CufflinksValue(String chr, int start, int end) {
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

    @Override
    public String getChr() {
        return chr;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getStart() {
        return start;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getEnd() {
        return end;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setStart(int start) {
        this.start = start;
    }

    @Override
    public void setEnd(int end) {
        this.end = end;
    }
}
