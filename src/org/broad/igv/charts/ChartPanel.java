package org.broad.igv.charts;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;
import java.io.Serializable;

/**
 * @author Jim Robinson
 * @date 10/26/11
 */
public class ChartPanel extends JPanel implements Serializable {

    ScatterPlot scatterPlot;

    private static final float[][] dash = {null, {1.0f, 1.0f}, {3.0f, 1.0f}, {4.0f, 4.0f}, {4.0f, 4.0f, 2.0f, 4.0f}};
    public static final BasicStroke DOT1 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, dash[1], 0.0f);
    public static final BasicStroke DOT2 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, dash[2], 0.0f);


    public ChartPanel() {
        init();
    }

    public void init() {
        this.setLayout(new BorderLayout());

        PlotPanel plotPanel = new PlotPanel();
        plotPanel.setBackground(Color.lightGray);
        plotPanel.setBorder(BorderFactory.createLineBorder(Color.BLACK, 1));
        this.add(plotPanel, BorderLayout.CENTER);

        AxisPanel xAxisPanel = new AxisPanel(AxisPanel.Orientation.HORIZONTAL);
        xAxisPanel.setPreferredSize(new Dimension(1000, 20));
        xAxisPanel.setBackground(Color.BLUE);
        this.add(xAxisPanel, BorderLayout.SOUTH);

        AxisPanel yAxisPanel = new AxisPanel(AxisPanel.Orientation.VERTICAL);
        yAxisPanel.setPreferredSize(new Dimension(20, 1000));
        yAxisPanel.setBackground(Color.GREEN);
        add(yAxisPanel, BorderLayout.WEST);

        LegendPanel legendPanel = new LegendPanel();
        legendPanel.setPreferredSize(new Dimension(100, 1000));
        legendPanel.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        add(legendPanel, BorderLayout.EAST);
    }

    public void setScatterPlotModel(ScatterPlot scatterPlot) {
        this.scatterPlot = scatterPlot;
        repaint();
    }

    class PlotPanel extends JPanel {

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            if (scatterPlot != null) {

                Rectangle r = new Rectangle(0, 0, getWidth(), getHeight());

                //       getVisibleRect();
                // todo -- use damager rectangle
                scatterPlot.draw((Graphics2D) g, r);
            }
        }
    }

    static class AxisPanel extends JComponent {
        enum Orientation {HORIZONTAL, VERTICAL}

        ;

        Orientation orientation;

        Axis axis;

        AxisPanel(Orientation orientation) {
            this.orientation = orientation;

        }

        public void setAxisModel(Axis axis) {
            this.axis = axis;
        }

        @Override
        protected void paintComponent(Graphics g) {
            g.setColor(getBackground());
            g.fillRect(0, 0, getWidth(), getHeight());
        }
    }

    static class LegendPanel extends JComponent {

    }
}
