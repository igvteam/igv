/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.data;

import org.igv.feature.LocusScore;
import org.igv.track.WindowFunction;

/**
 * @author jrobinso
 */
public class BasicScore implements LocusScore {

    protected int start;
    protected int end;
    protected float score;

    public BasicScore(int start, int end, float score) {

        this.start = start;
        this.end = end;
        this.score = score;
    }

    public BasicScore(BasicScore bs) {
        this.start = bs.start;
        this.end = bs.end;
        this.score = bs.score;
    }

    public BasicScore copy() {
        return new BasicScore(this);
    }

    /**
     * This method is required by the Tribble interface but not used.  To save space the chromosome is not stored, return null
     * @return
     */
    public String getChr() {
        return null;
    }

    @Override
    public String getContig() {
        return null;
    }

    public int getStart() {
        return start;
    }

    public float getScore() {
        return score;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    @Override
    public String toString() {
        return String.format("BasicScore: %d-%d ; %f", getStart(), getEnd(), getScore());
    }

    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        StringBuffer buf = new StringBuffer();
        buf.append(String.format("Value: %g at position %d",  score, (int)position));
        if(windowFunction != null) {
            buf.append("<br>Window function: " + windowFunction);
        }
        return buf.toString();
    }

}
