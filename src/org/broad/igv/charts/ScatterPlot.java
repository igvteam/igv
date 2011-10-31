package org.broad.igv.charts;

import org.broad.igv.track.AttributeManager;

import java.awt.*;
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
    Axis xAxis = new Axis(Axis.Orientation.HORIZONTAL);
    Axis yAxis = new Axis(Axis.Orientation.VERTICAL);
    XYDataModel dataModel;

    Rectangle pointShape = new Rectangle(7, 7);
    int offsetX = pointShape.getBounds().width / 2;
    int offsetY = pointShape.getBounds().height / 2;

    public void setModel(XYDataModel dataModel) {
        this.dataModel = dataModel;

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
        xAxis.setLabel(dataModel.getxLabel());
        yAxis.setLabel(dataModel.getyLabel());
    }

    public void draw(Graphics2D graphics, Rectangle bounds) {

        if (dataModel == null) return;

        // graphics.setColor(Color.white);
        // X ticks
        Color c = graphics.getColor();
        Stroke s = graphics.getStroke();

        drawGrid(graphics, bounds);

        graphics.setColor(c);
        graphics.setStroke(s);

        String categoryName = dataModel.getCategoryName();

        for (String sn : dataModel.getSeriesNames()) {

            Color color = getColor(categoryName, sn);

            graphics.setColor(color);

            List<XYDataPoint> dataPoints = dataModel.getDataPoints(sn);

            for (XYDataPoint dataPoint : dataPoints) {
                double x = dataPoint.getX();
                double y = dataPoint.getY();
                if (!Double.isNaN(x) && !Double.isNaN(y)) {

                    int px = xAxis.getPixelForValue(x);
                    final int bottom = bounds.y + bounds.height;
                    int pY = yAxis.getPixelForValue(y);
                    if (px >= bounds.x && px <= bounds.x + bounds.width &&
                            pY >= bounds.y && pY <= bottom) {
                        graphics.fillRect(px - offsetX, pY - offsetY, pointShape.width, pointShape.height);
                    }
                }
            }
        }
    }


    public XYDataPoint getDataPointAtPixel(int px, int py) {

        double x = xAxis.getDataValueForPixel(px);
        double y = yAxis.getDataValueForPixel(py);
        double toleranceX = ((pointShape.width + 1) / xAxis.getScale()) / 2;
        double toleranceY = ((pointShape.height + 1) / yAxis.getScale()) / 2;
        return dataModel.getDataPointAtPixel(x, y, toleranceX, toleranceY);
    }


    public static Color getColor(String categoryName, String sn) {
        Color color = (categoryName == null || categoryName.equals("")) ? Color.blue :
                AttributeManager.getInstance().getColor(categoryName, sn);

        // White is the "no-value" color in the attribute panel, but it doesn't work well on the plot. Switch to black
        if (color == Color.white) color = Color.black;
        return color;
    }

    private void drawGrid(Graphics2D graphics, Rectangle bounds) {

        graphics.setColor(VERY_LIGHT_GRAY);
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
        if (px > bounds.x && px < bounds.x + bounds.width) {
            graphics.drawLine(px, bounds.y, px, bounds.y + bounds.height);
        }
        py = yAxis.getPixelForValue(0);
        if (py > bounds.y && py < bounds.y + bounds.height) {
            graphics.drawLine(bounds.x, py, bounds.x + bounds.width, py);
        }


    }


}
