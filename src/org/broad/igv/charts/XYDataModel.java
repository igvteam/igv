package org.broad.igv.charts;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.util.ColorUtilities;

import java.awt.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 10/25/11
 */
public class XYDataModel {


    Map<String, XYSeries> seriesMap = new HashMap();
    String categoryName;
    String xLabel;
    String yLabel;

    public XYDataModel(String categoryName, String xLabel, String yLabel) {
        this.categoryName = categoryName;
        this.xLabel = xLabel;
        this.yLabel = yLabel;
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

    public String getxLabel() {
        return xLabel;
    }

    public String getyLabel() {
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
}
