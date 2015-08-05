/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;

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

    public String getValueString(double position, WindowFunction windowFunction) {
        StringBuffer buf = new StringBuffer();
        buf.append(String.format("Value: %g at position %d",  score, (int)position));
        if(windowFunction != null) {
            buf.append("<br>Window function: " + windowFunction);
        }
        return buf.toString();
    }

}
