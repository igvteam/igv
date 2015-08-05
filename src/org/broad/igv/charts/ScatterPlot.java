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

import org.broad.igv.PreferenceManager;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.awt.geom.Path2D;
import java.util.*;
import java.util.List;

/**
 * Renders a scatterplot
 *
 * @author Jim Robinson
 * @date 10/25/11
 */
public class ScatterPlot {


    public static final Color VERY_LIGHT_GRAY = new Color(250, 250, 250);

    PaletteColorTable colorTable;

    ScatterPlotData spData;
    Axis xAxis = new Axis(Axis.Orientation.HORIZONTAL);
    Axis yAxis = new Axis(Axis.Orientation.VERTICAL);
    private XYDataModel dataModel;
    Set<XYDataPoint> selectedPoints;

    Rectangle pointShape = new Rectangle(7, 7);
    int offsetX = pointShape.getBounds().width / 2;
    int offsetY = pointShape.getBounds().height / 2;

    HashSet<String> filteredSeries = new HashSet<String>();

    public static boolean isDataCategory(String selectedCategory) {
        return selectedCategory.equals(TrackType.COPY_NUMBER.toString()) ||
                selectedCategory.equals(TrackType.GENE_EXPRESSION.toString()) ||
                selectedCategory.equals(TrackType.DNA_METHYLATION.toString());
    }


    public ScatterPlot(ScatterPlotData spData) {
        this.spData = spData;

    }

    public synchronized void setModel(XYDataModel dataModel) {

        this.dataModel = dataModel;

        colorTable = new PaletteColorTable(ColorUtilities.getPalette("Set 1"));
        double minX = Double.MAX_VALUE;
        double maxX = -minX;
        double minY = minX;
        double maxY = maxX;

        for (String sn : dataModel.getSeriesNames()) {
            List<XYDataPoint> dataPoints = dataModel.getDataPoints(sn);
            for (XYDataPoint dataPoint : dataPoints) {
                double x = dataPoint.getX();
                double y = dataPoint.getY();
                if (!Double.isNaN(x) && !Double.isNaN(y)) {
                    minX = Math.min(minX, x);
                    maxX = Math.max(maxX, x);
                    minY = Math.min(minY, y);
                    maxY = Math.max(maxY, y);
                }
            }
        }
        xAxis.setRange(minX, maxX);
        yAxis.setRange(minY, maxY);
        xAxis.setLabel(dataModel.getXLabel());
        yAxis.setLabel(dataModel.getYLabel());
    }

    public void draw(Graphics2D graphics, Rectangle bounds, Rectangle clipRect) {

        if (dataModel == null) return;

        // graphics.setColor(Color.white);
        // X ticks
        Color c = graphics.getColor();
        Stroke s = graphics.getStroke();

        drawGrid(graphics, bounds, clipRect);

        graphics.setColor(c);
        graphics.setStroke(s);

        String categoryName = dataModel.getCategoryName();

        for (String sn : dataModel.getSeriesNames()) {

            if (filteredSeries.contains(sn)) continue;

            Color color = null;
            double[] categoryValues = null;

            if (isDataCategory(categoryName)) {
                categoryValues = spData.getDataValues(categoryName);
            } else {
                color = getColor(categoryName, sn);
                graphics.setColor(color);
            }


            List<XYDataPoint> dataPoints = dataModel.getDataPoints(sn);

            for (XYDataPoint dataPoint : dataPoints) {
                double x = dataPoint.getX();
                double y = dataPoint.getY();
                if (!Double.isNaN(x) && !Double.isNaN(y)) {
                    int px = xAxis.getPixelForValue(x);
                    int pY = yAxis.getPixelForValue(y);
                    if (clipRect.contains(px, pY)) {

                        if (isDataCategory(categoryName)) {
                            color = getDataColor(categoryName, categoryValues[dataPoint.getIdx()]);
                            graphics.setColor(color);
                        }

                        graphics.fillOval(px - offsetX, pY - offsetY, pointShape.width, pointShape.height);
                    }
                }
            }

            // Mutations
            Stroke stroke = graphics.getStroke();
            BasicStroke thickLine = new BasicStroke(2);
            graphics.setStroke(thickLine);
            graphics.setColor(Color.black);
            for (XYDataPoint dataPoint : dataPoints) {
                double x = dataPoint.getX();
                double y = dataPoint.getY();
                if (!Double.isNaN(x) && !Double.isNaN(y)) {
                    int px = xAxis.getPixelForValue(x);
                    int pY = yAxis.getPixelForValue(y);
                    if (clipRect.contains(px, pY) && dataPoint.getMutationCount() > 0) {
                        graphics.drawOval(px - offsetX - 1, pY - offsetY - 1, pointShape.width + 1, pointShape.height + 1);
                    }
                }
            }
            graphics.setStroke(stroke);
            // Outline selected points.  this is done here, inside the "series" loop, so that series filtering is
            // respected.  I would be more effecient to do it outside
            if (selectedPoints != null) {
                graphics.setColor(Color.ORANGE);
                for (XYDataPoint dataPoint : dataPoints) {
                    if (selectedPoints.contains(dataPoint)) {
                        double x = dataPoint.getX();
                        double y = dataPoint.getY();
                        if (!Double.isNaN(x) && !Double.isNaN(y)) {
                            int px = xAxis.getPixelForValue(x);
                            int pY = yAxis.getPixelForValue(y);
                            if (clipRect.contains(px, pY)) {
                                graphics.drawRect(px - offsetX, pY - offsetY, pointShape.width, pointShape.height);
                            }
                        }
                    }
                }
            }
        }

        // Outline selected points
//        graphics.setColor(Color.ORANGE);
//        if (selectedPoints != null) {
//            for (XYDataPoint dataPoint : selectedPoints) {
//                double x = dataPoint.getX();
//                double y = dataPoint.getY();
//                if (!Double.isNaN(x) && !Double.isNaN(y)) {
//                    int px = xAxis.getPixelForValue(x);
//                    int pY = yAxis.getPixelForValue(y);
//                    if (clipRect.contains(px, pY)) {
//                        graphics.drawRect(px - offsetX - 1, pY - offsetY - 1, pointShape.width + 2, pointShape.height + 2);
//                    }
//                }
//
//            }
//        }
    }


    Map<TrackType, ContinuousColorScale> colorScales = new HashMap<TrackType, ContinuousColorScale>();

    private Color getDataColor(String categoryName, double categoryValue) {
        TrackType tt = TrackType.valueOf(categoryName);
        ContinuousColorScale scale = colorScales.get(tt);
        if (scale == null) {
            scale = PreferenceManager.getInstance().getColorScale(tt);// IGV.getInstance().getSession().getColorScale(tt);
        }
        return scale.getColor((float) categoryValue);

    }

    public XYDataPoint getDataPointAtPixel(int px, int py) {

        double x = xAxis.getDataValueForPixel(px);
        double y = yAxis.getDataValueForPixel(py);
        double toleranceX = ((pointShape.width + 1) / xAxis.getScale()) / 2;
        double toleranceY = ((pointShape.height + 1) / yAxis.getScale()) / 2;
        return dataModel.getDataPointAtPixel(x, y, toleranceX, toleranceY);
    }


    public Color getColor(String categoryName, String sn) {

        return colorTable.get(sn);

//        Color color;
//        if (categoryName == null || categoryName.equals("")) {
//            color = Color.blue;
//        } else if (categoryName.equals("Mut Count")) {
//            if (sn == null || sn.equals("") || sn.equals("Unknown")) {
//                color = Color.darkGray;
//            } else if (sn.equals("0")) {
//                color = Color.green.darker();
//            } else if (sn.equals("1")) {
//                color = Color.blue;
//            } else if (sn.equals("2")) {
//                color = Color.orange;
//            } else {
//                color = Color.red;
//            }
//        } else {
//            color = AttributeManager.getInstance().getColor(categoryName, sn);
//        }
//
//        // White is the "no-value" color in the attribute panel, but it doesn't work well on the plot. Switch to black
//        if (color == Color.white) color = Color.darkGray;
//        return color;
    }

    private void drawGrid(Graphics2D graphics, Rectangle bounds, Rectangle clipRect) {

        graphics.setColor(new Color(161, 196, 214));
        graphics.setStroke(ChartPanel.DOT1);
        double[] xticks = xAxis.ticks;
        double xtick = xticks[0];
        int px = 0;
        while (px < bounds.x + bounds.width) {
            px = xAxis.getPixelForValue(xtick);
            if (px > bounds.x && px < bounds.x + bounds.width) {
                graphics.drawLine(px, bounds.y, px, bounds.y + bounds.height);
            }
            xtick += xticks[1];
        }
        double[] yticks = yAxis.ticks;
        double ytick = yticks[0];
        int py = bounds.y + bounds.height;
        while (py > bounds.y) {
            py = yAxis.getPixelForValue(ytick);
            if (py > bounds.y && py < bounds.y + bounds.height) {
                graphics.drawLine(bounds.x, py, bounds.x + bounds.width, py);
            }
            ytick += yticks[1];
        }

        // Emphasize zero
        graphics.setColor(Color.blue.darker());
        graphics.setStroke(ChartPanel.DOT2);
        px = xAxis.getPixelForValue(0);
        if (px > clipRect.x && px < clipRect.x + clipRect.width) {
            graphics.drawLine(px, clipRect.y, px, clipRect.y + clipRect.height);
        }
        py = yAxis.getPixelForValue(0);
        if (py > clipRect.y && py < clipRect.y + clipRect.height) {
            graphics.drawLine(clipRect.x, py, clipRect.x + clipRect.width, py);
        }


    }


    public void selectPointsInPath(Path2D path) {
        if (dataModel != null) {
            selectedPoints = dataModel.getDataPointsIn(path);
        }
    }

    public void clearSelections() {
        selectedPoints = null;
    }

    public void addSeriesFilter(String sn) {
        filteredSeries.add(sn);
    }

    public void removeSeriesFilter(String sn) {
        filteredSeries.remove(sn);
    }

    public XYDataModel getDataModel() {
        return dataModel;
    }
}
