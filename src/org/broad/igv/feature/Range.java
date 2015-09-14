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

package org.broad.igv.feature;

import htsjdk.tribble.Feature;
import org.broad.igv.util.Utilities;


/**
 * Basic class to specify a genomic interval.
 * Coordinates are intended to be 0-based half-open
 * Please do not add or remove any fields (want to keep it very simple).
 * Additional methods for calculating overlap are okay
 *
 * @author jacob
 * @author 2013-May-20
 */
public class Range implements Feature {

    protected String chr = null;
    protected int start = -1;
    protected int end = -1;

    public Range(String chr, int start, int end){
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

    public int getLength(){
        return end - start;
    }

    /**
     * Determine whether this interval fully contains the specified
     * input interval.
     *
     * A negative input start position has special meaning.  It is considered within the interval if the interval
     * contains position "0".
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public boolean contains(String chr, int start, int end) {
        return Utilities.objectEqual(this.chr, chr)
                && this.start <= (start < 0 ? 0 : start)
                && this.end >= end;
    }

    /**
     * Determine whether there is any overlap between this interval and the specified interval
     *
     * Negative positions have no special meaning
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public boolean overlaps(String chr, int start, int end) {
        return Utilities.objectEqual(this.chr, chr) && this.start <= end && this.end >= start;
    }

    public boolean overlaps(Range range) {
        return this.overlaps(range.getChr(), range.getStart(), range.getEnd());
    }

}
