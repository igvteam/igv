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

import java.awt.geom.Path2D;
import java.util.*;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 10/25/11
 */
public class XYDataModel {


   private Map<String, XYSeries> seriesMap = new HashMap();
    private String categoryName;
    private String xLabel;
    private String yLabel;
    ScatterPlotData scatterPlotData;


    public XYDataModel(String categoryName, String xLabel, String yLabel, ScatterPlotData scatterPlotData) {
        this.categoryName = categoryName;
        this.xLabel = xLabel;
        this.yLabel = yLabel;
        this.scatterPlotData = scatterPlotData;
    }


    public void addSeries(XYSeries xySeries) {

        int idx = seriesMap.size();
        seriesMap.put(xySeries.getName(), xySeries);
    }

    public int numPoints(String seriesName) {
        XYSeries series = seriesMap.get(seriesName);
        return series == null ? 0 : series.getSize();

    }

    public List<XYDataPoint> getDataPoints(String seriesName) {
        XYSeries series = seriesMap.get(seriesName);
        return series == null ? null : series.getDataPoints();
    }

    public Collection<String> getSeriesNames() {
        return seriesMap.keySet();
    }

    public String getCategoryName() {
        return categoryName;
    }

    public String getXLabel() {
        return xLabel;
    }

    public String getYLabel() {
        return yLabel;
    }

    /**
     * Return the first data point found that contains the given point.
     * TODO -- impose z-order, last-one-drawn, or some other tie-breaker?
     *
     * @param x
     * @param y
     * @param toleranceX
     * @param toleranceY
     * @return
     */
    public XYDataPoint getDataPointAtPixel(double x, double y, double toleranceX, double toleranceY) {

        for (XYSeries series : seriesMap.values()) {
            XYDataPoint dp = series.getDataPoint(x, y, toleranceX, toleranceY);
            if (dp != null) {
                return dp;
            }
        }

        return null;
    }

    public Set<XYDataPoint> getDataPointsIn(Path2D path) {

        HashSet<XYDataPoint> points = new HashSet<XYDataPoint>();
        for (XYSeries series : seriesMap.values()) {
            Collection<XYDataPoint> dataPoints = series.getDataPoints();
            for(XYDataPoint dataPoint : dataPoints) {
                if(path.contains(dataPoint.getX(), dataPoint.getY())) {
                    points.add(dataPoint);
                }
            }
        }
        return points;


    }

}
