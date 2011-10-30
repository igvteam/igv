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

    public double[] getX(String seriesName) {
        XYSeries series = seriesMap.get(seriesName);
        return series == null ? null : series.getX();
    }

    public double[] getY(String seriesName) {
        XYSeries series = seriesMap.get(seriesName);
        return series == null ? null : series.getY();
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
}
