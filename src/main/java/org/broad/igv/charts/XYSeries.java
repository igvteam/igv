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


