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

package org.broad.igv.feature.exome;


import htsjdk.tribble.Feature;

/**
 * @author Jim Robinson
 * @date 5/24/12
 */
public class ExomeBlock implements Feature {
    int idx;
    private int genomeStart;
    private int exomeStart;
    private int length;
    private int leftPixel;
    private int rightPixel;

    public ExomeBlock(int idx, int genomeStart, int xomeStart, int length) {
        this.idx = idx;
        this.genomeStart = genomeStart;
        this.exomeStart = xomeStart;
        this.length = length;
    }

    public void setScreenBounds(int leftPixel, int rightPixel) {
        this.leftPixel = leftPixel;
        this.rightPixel = rightPixel;
    }

    public int getLeftPixel() {
        return leftPixel;
    }

    public int getRightPixel() {
        return rightPixel;
    }

    public int getIdx() {
        return idx;
    }

    public int getGenomeStart() {
        return genomeStart;
    }

    public int getGenomeEnd() {
        return genomeStart + length;
    }

    public int getExomeStart() {
        return exomeStart;
    }

    public int getExomeEnd() {
        return exomeStart + length;
    }

    public void extend(int x) {
        if (x > genomeStart + length) {
            length = x - genomeStart;
        };
    }

    public int getLength() {
        return getGenomeEnd() - getGenomeStart();
    }


    public String toString() {
        return "Block " + idx + " [" + genomeStart + ", " + getGenomeEnd() + ", " + exomeStart + ", " + length +"]";
    }

    @Override
    public String getChr() {

        return null;
    }

    @Override
    public String getContig() {
        return null;
    }

    @Override
    public int getStart() {
        return genomeStart;
    }

    @Override
    public int getEnd() {
        return genomeStart + length;
    }

    public int compareGenomePosition(double genomicPosition) {
        if(genomicPosition < genomeStart) {
            return -1;
        }
        else if(genomicPosition >= genomeStart + length) {
            return 1;
        }
        else {
            return 0;
        }

    }

    public int compareExomePosition(int exomePosition) {
        if(exomePosition < exomeStart) {
            return -1;
        }
        else if(exomePosition >= (exomeStart + length)) {
            return 1;
        }
        else {
            return 0;
        }
    }
}
