package org.broad.igv.charts;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.FontManager;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
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

        xAxisPanel = new AxisPanel();
        xAxisPanel.setPreferredSize(new Dimension(1000, 60));
        xAxisPanel.setFont(FontManager.getDefaultFont());
        this.add(xAxisPanel, ChartLayout.XAXIS);


        yAxisPanel = new AxisPanel();
        yAxisPanel.setPreferredSize(new Dimension(60, 1000));
        yAxisPanel.setFont(FontManager.getDefaultFont());
        add(yAxisPanel, ChartLayout.YAXIS);

        legendPanel = new LegendPanel();
        legendPanel.setPreferredSize(new Dimension(150, 1000));
        legendPanel.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        add(legendPanel, ChartLayout.LEGEND);
    }

    public void setScatterPlotModel(ScatterPlot scatterPlot) {
        this.scatterPlot = scatterPlot;
        xAxisPanel.setAxisModel(scatterPlot.xAxis);
        yAxisPanel.setAxisModel(scatterPlot.yAxis);
        repaint();
    }

    class PlotPanel extends JPanel {

        PlotPanel() {

            setToolTipText("Plot panel");

            this.addMouseMotionListener(new MouseAdapter() {
                @Override
                public void mouseMoved(MouseEvent mouseEvent) {
                    if (scatterPlot != null) {
                        XYDataPoint dp = scatterPlot.getDataPointAtPixel(mouseEvent.getX(), mouseEvent.getY());
                        if (dp != null) {
                            String description = dp.getDescription();
                            if (description != null) {
                                PlotPanel.this.setToolTipText(description);
                                System.out.println(description);
                            }
                        } else {
                            PlotPanel.this.setToolTipText("");
                        }

                    }
                }
            });
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            if (scatterPlot != null) {

                Rectangle r = new Rectangle(0, 0, getWidth(), getHeight());

                //       getVisibleRect();
                // todo -- use damage rectangle
                scatterPlot.draw((Graphics2D) g, r);
            }
        }


    }

    class LegendPanel extends JComponent {

        int leftMargin = 10;
        int topMargin = 30;
        private Font labelFont;
        private Font headerFont;

        LegendPanel() {
            Font defaultFont = FontManager.getDefaultFont();
            labelFont = defaultFont.deriveFont(12);
            headerFont = defaultFont.deriveFont(14);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            Graphics2D g2D = (Graphics2D) g;

            if (scatterPlot != null) {
                Color color = g.getColor();
                Font font = g.getFont();

                Rectangle pointShape = scatterPlot.pointShape;

                String categoryName = scatterPlot.dataModel.categoryName;
                if (categoryName == null || categoryName.equals("")) {
                    return;
                }

                g2D.setFont(headerFont);
                g2D.drawString(categoryName, leftMargin, topMargin);

                g2D.setFont(labelFont);
                int y = topMargin + 20;
                for (String sn : scatterPlot.dataModel.getSeriesNames()) {
                    Color c = ScatterPlot.getColor(categoryName, sn);
                    g2D.setColor(c);
                    g2D.fillRect(leftMargin + 5, y, pointShape.width, pointShape.height);

                    String displayString = sn.equals("") ? "Unknown" : sn;
                    g2D.setColor(Color.black);
                    g2D.drawString(displayString, leftMargin + 20, y + 5);

                    y += 20;

                }


                g2D.setColor(color);
                g2D.setFont(font);


            }


        }

    }
}
