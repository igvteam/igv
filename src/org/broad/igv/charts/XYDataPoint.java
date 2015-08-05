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

package org.broad.igv.charts;

/**
 * @author Jim Robinson
 * @date 10/31/11
 */
public class XYDataPoint {

    private int idx;
    private double x;
    private double y;
    private int mutationCount;   // <= TODO make this some generic attribute, we've lost the generality of this class
    private String description;

    public XYDataPoint(int idx, double x, double y, int mutationCount, String description) {
        this.idx = idx;
        this.x = x;
        this.y = y;
        this.mutationCount = mutationCount;
        this.description = description;
    }


    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public int getMutationCount() {
        return mutationCount;
    }

    public String getDescription() {
        return description;
    }

    /**
     * Test if the point represented by (px, py) is with tolerance of this point.
     *
     *
     * @param px
     * @param py
     * @param toleranceX
     * @param toleranceY
     * @return
     */
    public boolean contains(double px, double py, double toleranceX, double toleranceY) {

        return px > x - toleranceX && px < x + toleranceX && py > y - toleranceY && py < y + toleranceY;

    }

    public int getIdx() {
        return idx;
    }
}
