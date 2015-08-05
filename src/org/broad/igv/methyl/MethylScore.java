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

package org.broad.igv.methyl;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.WindowFunction;

/**
 * @author Jim Robinson
 * @date 4/19/12
 */
public class MethylScore implements LocusScore {

    String chr;
    int start;
    int end;
    Strand strand;
    float percentMethylated;
    int totalCount;

    public MethylScore(String chr, int start, int end, Strand strand, float percentMethylated, int totalCount) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.percentMethylated = percentMethylated;
        this.totalCount = totalCount;
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

    public Strand getStrand() {
        return strand;
    }

    public float getScore() {
        return percentMethylated;
    }

    public int getCount() {
        return totalCount;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    /**
     * Return a string to be used for popup text.   The WindowFunction is passed
     * in so it can be used t annotate the value.  The LocusScore object itself
     * does not "know" from what window function it was derived
     *
     * @param position       Zero-based genome position
     * @param windowFunction
     * @return
     */
    public String getValueString(double position, WindowFunction windowFunction) {
        return percentMethylated + "%" + " [" + totalCount + "]" +
                (strand == Strand.POSITIVE ? " (+)" : (strand == Strand.NEGATIVE ? " (-)" : ""));
    }

}
