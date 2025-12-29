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
