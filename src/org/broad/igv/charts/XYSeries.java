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

import org.broad.igv.util.collections.DoubleArrayList;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 10/26/11
 */
public class XYSeries {

    String name;
    List<XYDataPoint> dataPoints;

    public XYSeries(String name) {
        this.name = name;
        dataPoints = new ArrayList<XYDataPoint>(500);
    }

    public void add(int idx, double x, double y, int mutationCount, String description) {
        dataPoints.add(new XYDataPoint(idx, x, y, mutationCount, description));
    }

    public int getSize() {
        return dataPoints.size();
    }

    public String getName() {
        return name;
    }

    public List<XYDataPoint> getDataPoints() {
        return dataPoints;
    }

    public XYDataPoint getDataPoint(double x, double y, double toleranceX, double toleranceY) {

        for (XYDataPoint dp : dataPoints) {
            if (dp.contains(x, y, toleranceX, toleranceY)) {
                return dp;
            }
        }
        return null;
    }
}


