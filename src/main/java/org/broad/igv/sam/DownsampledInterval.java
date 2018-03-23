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

package org.broad.igv.sam;

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
