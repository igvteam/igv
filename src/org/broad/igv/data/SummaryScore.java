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
public class SummaryScore implements LocusScore {

    int start;
    int end;
    float value;

    public SummaryScore(int start, int end, float value) {
        this.start = start;
        this.end = end;
        this.value = value;
    }

    public SummaryScore(LocusScore anotherScore) {
        this.start = anotherScore.getStart();
        this.end = anotherScore.getEnd();
        this.value = anotherScore.getScore();
    }

    public SummaryScore copy() {
        return new SummaryScore(this);
    }

    public void setStart(int start) {
        this.start = start;
    }

    public String getChr() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
    @Override
    public String getContig() {
        return null;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public float getValue() {
        return value;
    }

    public void setScore(float score) {
        this.value = score;
    }

    public float getScore() {
        return value;
    }

    public String getValueString(double position, WindowFunction wf) {
        return "Value:  " + getScore() + (wf == null ? "" : " (" + wf.toString() + ")");
    }

}
