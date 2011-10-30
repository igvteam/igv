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

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat();
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

    static class AxisPanel extends JComponent {

        private static final int TICK_SIZE = 4;
        private static final int TICK_GAP = 2;

        enum Orientation {HORIZONTAL, VERTICAL}

        Orientation orientation;
        Axis axis;

        AxisPanel(Orientation orientation) {
            this.orientation = orientation;

        }

        public void setAxisModel(Axis axis) {
            this.axis = axis;
        }


        @Override
        public void setBounds(int x, int y, int w, int h) {
            super.setBounds(x, y, w, h);
            updateAxisDimension(w, h);
        }

        private void updateAxisDimension(int w, int h) {

            if (orientation == Orientation.HORIZONTAL) {
                 axis.setPanelSize(w);
             } else {
                 axis.setPanelSize(h);
             }
         }

        @Override
        protected void paintComponent(Graphics g) {


            final FontMetrics fontMetrics = g.getFontMetrics();
            int minusStrWidth = fontMetrics.stringWidth("-") / 2;
            int strHeight = fontMetrics.getHeight();
            int bottom = getHeight();

            if (orientation == Orientation.HORIZONTAL) {
                double[] xticks = axis.ticks;
                double xtick = xticks[0];
                int px = 0;
                final int width = getWidth();
                final int tickLabelY = TICK_SIZE + 2 * TICK_GAP + strHeight;
                while (px < width) {
                    px = axis.getPixelForValue(xtick);
                    if (px > 0 && px < width) {
                        g.drawLine(px, TICK_GAP, px, TICK_GAP + TICK_SIZE);

                        String label = DECIMAL_FORMAT.format(xtick);
                        int strWidth = fontMetrics.stringWidth(label);
                        if (px > strWidth && (px + strWidth) < width) {
                            int strPosition = px - strWidth / 2;
                            if (xtick < 0) {
                                strPosition -= minusStrWidth;
                            }
                            g.drawString(label, strPosition, tickLabelY);
                        }

                    }
                    xtick += xticks[1];
                }

                // If there is room draw the label
                if (bottom - strHeight - 5 > tickLabelY) {
                    String label = axis.getLabel();
                    if (label != null) {
                        int strWidth = fontMetrics.stringWidth(label);
                        int strX = (getWidth() - strWidth) / 2;
                        g.drawString(label, strX, bottom - 5);
                    }
                }


            }
        }
    }

    static class LegendPanel extends JComponent {

    }
}
