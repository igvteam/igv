package org.broad.igv.charts;

import org.broad.igv.ui.FontManager;

import javax.swing.*;
import java.awt.*;
import java.io.Serializable;
import java.text.DecimalFormat;

/**
 * @author Jim Robinson
 * @date 10/26/11
 */
public class ChartPanel extends JPanel implements Serializable {

    ScatterPlotRenderer scatterPlotRenderer;

    private static final float[][] dash = {null, {1.0f, 1.0f}, {3.0f, 1.0f}, {4.0f, 4.0f}, {4.0f, 4.0f, 2.0f, 4.0f}};
    public static final BasicStroke DOT1 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, dash[1], 0.0f);
    public static final BasicStroke DOT2 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, dash[2], 0.0f);

    PlotPanel plotPanel;
    AxisPanel xAxisPanel;
    AxisPanel yAxisPanel;
    LegendPanel legendPanel;

    public ChartPanel() {
        init();
    }

    public void init() {
        this.setLayout(new ChartLayout());

        plotPanel = new PlotPanel();
        plotPanel.setBackground(Color.lightGray);
        plotPanel.setBorder(BorderFactory.createLineBorder(Color.BLACK, 1));
        this.add(plotPanel, ChartLayout.CHART);

        xAxisPanel = new AxisPanel(AxisPanel.Orientation.HORIZONTAL);
        xAxisPanel.setPreferredSize(new Dimension(1000, 60));
        xAxisPanel.setFont(FontManager.getDefaultFont());
        this.add(xAxisPanel, ChartLayout.XAXIS);

        xAxisPanel.setBorder(BorderFactory.createLineBorder(Color.green));

        yAxisPanel = new AxisPanel(AxisPanel.Orientation.VERTICAL);
        yAxisPanel.setPreferredSize(new Dimension(60, 1000));
        yAxisPanel.setFont(FontManager.getDefaultFont());
        add(yAxisPanel, ChartLayout.YAXIS);

        yAxisPanel.setBorder(BorderFactory.createLineBorder(Color.cyan));

        legendPanel = new LegendPanel();
        legendPanel.setPreferredSize(new Dimension(100, 1000));
        legendPanel.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        add(legendPanel, ChartLayout.LEGEND);
    }

    public void setScatterPlotModel(ScatterPlotRenderer scatterPlotRenderer) {
        this.scatterPlotRenderer = scatterPlotRenderer;
        xAxisPanel.setAxisModel(scatterPlotRenderer.xAxis);
        yAxisPanel.setAxisModel(scatterPlotRenderer.yAxis);
        repaint();
    }

    class PlotPanel extends JPanel {

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            if (scatterPlotRenderer != null) {

                Rectangle r = new Rectangle(0, 0, getWidth(), getHeight());

                //       getVisibleRect();
                // todo -- use damager rectangle
                scatterPlotRenderer.draw((Graphics2D) g, r);
            }
        }
    }

    static class LegendPanel extends JComponent {

    }
}
