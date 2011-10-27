package org.broad.igv.charts;

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

    static Color [] color = {Color.blue, Color.green, Color.yellow};

    Map<String, XYSeries> seriesMap = new HashMap();
    Map<String, Color> seriesColorMap = new HashMap<String, Color>();

    public XYDataModel() {
        //To change body of created methods use File | Settings | File Templates.
    }

    public void addSeries(XYSeries xySeries) {

        int idx = seriesMap.size();
        Color c = idx >= color.length ? ColorUtilities.randomColor(idx) : color[idx];
        seriesColorMap.put(xySeries.getName(), c);

        seriesMap.put(xySeries.getName(), xySeries);
    }

    public int numPoints(String seriesName) {
        XYSeries series =  seriesMap.get(seriesName);
        return series == null ? 0 : series.getSize();

    }
    public double[] getX(String seriesName) {
        XYSeries series =  seriesMap.get(seriesName);
        return series == null ? null : series.getX();
    }

    public double[] getY(String seriesName) {
        XYSeries series =  seriesMap.get(seriesName);
        return series == null ? null : series.getY();
    }

    public Color getColor(String seriesName) {
        return seriesColorMap.get(seriesName);
    }

    public Collection<String> getSeriesNames() {
        return seriesMap.keySet();
    }
}
