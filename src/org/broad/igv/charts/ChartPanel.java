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
import org.broad.igv.ui.FontManager;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Path2D;
import java.io.Serializable;
import java.util.*;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 10/26/11
 */
public class ChartPanel extends JPanel implements Serializable {

    ScatterPlot scatterPlot;

    public static final BasicStroke DOT1 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER,
            10.0f, new float[]{1.0f, 1.0f}, 0.0f);
    public static final BasicStroke DOT2 = new BasicStroke(1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER,
            10.0f, new float[]{3.0f, 1.0f}, 0.0f);

    PlotPanel plotPanel;
    AxisPanel xAxisPanel;
    AxisPanel yAxisPanel;
    LegendPanel legendPanel;
//    ToolBar toolPanel;

    boolean lassoInProgress = false;
    SelectionPath lassoPath = null;

    public ChartPanel() {
        init();
    }

    public void init() {
        this.setLayout(new ChartLayout());


//        toolPanel = new ToolBar();
//        toolPanel.setPreferredSize(new Dimension(1000, 20));
//        this.add(toolPanel, ChartLayout.TITLE);


        plotPanel = new PlotPanel();
        plotPanel.setBackground(PreferenceManager.getInstance().getAsColor(PreferenceManager.BACKGROUND_COLOR));
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
        legendPanel.setPreferredSize(new Dimension(200, 1000));

        final Border lineBorder = BorderFactory.createLineBorder(Color.BLACK);

        legendPanel.setBorder(lineBorder);
        add(legendPanel, ChartLayout.LEGEND);
    }

    public void setScatterPlotModel(ScatterPlot scatterPlot) {
        this.scatterPlot = scatterPlot;
        xAxisPanel.setAxisModel(scatterPlot.xAxis);
        yAxisPanel.setAxisModel(scatterPlot.yAxis);
        legendPanel.rebuild();
        repaint();
    }

    class PlotPanel extends JPanel {

        PlotPanel() {

            setToolTipText("Plot panel");

            final MouseAdapter mouseAdapter = new MouseAdapter() {
                @Override
                public void mouseMoved(MouseEvent mouseEvent) {
                    if (scatterPlot != null) {
                        if (lassoInProgress) {
                            // Ignore
                        } else {
                            XYDataPoint dp = scatterPlot.getDataPointAtPixel(mouseEvent.getX(), mouseEvent.getY());
                            if (dp != null) {
                                String description = dp.getDescription();
                                if (description != null) {
                                    PlotPanel.this.setToolTipText(description);
                                }
                            } else {
                                PlotPanel.this.setToolTipText("");
                            }
                        }

                    }
                }

                @Override
                public void mouseDragged(MouseEvent mouseEvent) {
                    if (lassoInProgress) {
                        lassoPath.addPoint(mouseEvent.getPoint());
                        paintImmediately(new Rectangle(0, 0, getWidth(), getHeight())); //lassoPath.getBounds());
                    }
                }

                @Override
                public void mouseReleased(MouseEvent mouseEvent) {
                    if (lassoInProgress) {
                        lassoInProgress = false;
                        if (lassoPath.size() > 2) {

                            Axis xAxis = scatterPlot.xAxis;
                            Axis yAxis = scatterPlot.yAxis;

                            Path2D path = new Path2D.Double(java.awt.geom.Path2D.WIND_NON_ZERO, lassoPath.size());
                            Iterator<Point> iter = lassoPath.getPoints().iterator();
                            Point p = iter.next();
                            double x = xAxis.getDataValueForPixel(p.x);
                            double y = yAxis.getDataValueForPixel(p.y);
                            path.moveTo(x, y);
                            while (iter.hasNext()) {
                                p = iter.next();
                                x = xAxis.getDataValueForPixel(p.x);
                                y = yAxis.getDataValueForPixel(p.y);
                                path.lineTo(x, y);
                            }
                            path.closePath();
                            scatterPlot.selectPointsInPath(path);
                            Rectangle damageRect = lassoPath.getBounds();
                            lassoPath = null;
                            repaint(damageRect);


                        }

                    }
                }

                @Override
                public void mouseClicked(MouseEvent mouseEvent) {

                    if (lassoInProgress) {
                        lassoPath.addPoint(mouseEvent.getPoint());
                    } else {
                        // TODO -- usual click options, if unmodified and on a point select that, if modifier key
                        // TODO -- add the point, etc
                        scatterPlot.clearSelections();
                        repaint();
                    }

                }
            };

            this.addMouseListener(mouseAdapter);
            this.addMouseMotionListener(mouseAdapter);
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

             if (scatterPlot != null) {
                Rectangle r = new Rectangle(0, 0, getWidth(), getHeight());
                scatterPlot.draw((Graphics2D) g, r, g.getClipBounds());
            }

//            if (lassoPath != null && lassoPath.size() > 1) {
//                Iterator<Point> iter = lassoPath.getPoints().iterator();
//                Point lastPoint = iter.next();
//                g.setColor(Color.cyan);
//                while (iter.hasNext()) {
//                    Point point = iter.next();
//                    g.drawLine(lastPoint.x, lastPoint.y, point.x, point.y);
//                    lastPoint = point;
//                }
//            }
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
            headerFont = defaultFont.deriveFont(Font.BOLD, 16);

            final BoxLayout boxLayout = new BoxLayout(this, BoxLayout.Y_AXIS);

            setLayout(boxLayout);


            rebuild();

        }

        void rebuild() {
            if (scatterPlot == null || scatterPlot.getDataModel() == null) return;

            removeAll();


            String categoryName = scatterPlot.getDataModel().getCategoryName();
            if (categoryName == null || categoryName.equals("")) {
                return;
            }


            add(Box.createVerticalStrut(topMargin));

            JLabel catLabel = new JLabel(" " + categoryName);
            catLabel.setFont(headerFont);
            add(catLabel);
            add(new JLabel(" "));    // Spacer


            if (ScatterPlot.isDataCategory(categoryName)) {


            } else {
                Rectangle pointShape = new Rectangle(10, 10); //scatterPlot.pointShape;

                // Rebuild the "series" checkboxes
                java.util.List<String> seriesNames = new ArrayList<String>(scatterPlot.getDataModel().getSeriesNames());

                sortSeriesNames(seriesNames);

                for (final String sn : seriesNames) {
                    Color c = scatterPlot.getColor(categoryName, sn);
                    LegendIcon icon = new LegendIcon(pointShape, c);

                    String labelString = (sn.trim()).equals("") ? "<No Value>" : sn;

                    JLabel label = new JLabel(labelString, icon, SwingConstants.LEFT);
                    label.setFont(labelFont);

                    final JCheckBox cb = new JCheckBox();
                    cb.setSelected(true);
                    cb.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent actionEvent) {
                            if (cb.isSelected()) {
                                scatterPlot.removeSeriesFilter(sn);
                            } else {
                                scatterPlot.addSeriesFilter(sn);
                            }
                            plotPanel.repaint();
                        }
                    });

                    JPanel panel = new JPanel();
                    panel.setAlignmentX(LEFT_ALIGNMENT);
                    panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));

                    panel.add(cb);
                    panel.add(label);
                    add(panel);
                }
                revalidate();
            }

        }

        /**
         * Sort the list of series names so that
         * (1) numbers come before letters
         * (2) numbers are in ascending numeric order
         * (3) "blank" (empty string) comes last
         *
         * @param seriesNames
         */
        private void sortSeriesNames(List<String> seriesNames) {
            Collections.sort(seriesNames, new Comparator<String>() {
                public int compare(String s1, String s2) {
                    // Try numeric sort first
                    double d1;
                    try {
                        d1 = Double.parseDouble(s1);
                    } catch (NumberFormatException e) {
                        d1 = Double.MAX_VALUE;
                    }
                    double d2;
                    try {
                        d2 = Double.parseDouble(s2);
                    } catch (NumberFormatException e) {
                        d2 = Double.MAX_VALUE;
                    }
                    if (d1 == d2) {
                        // do standard comapre, but put blanks at end
                        if (s1.equals("")) s1 = "ZZZ";
                        if (s2.equals("")) s2 = "ZZZ";
                        return s1.compareTo(s2);
                    } else {
                        return (int) (d1 - d2);
                    }
                }
            });
        }

//        @Override
//        protected void paintComponent(Graphics g) {
//            super.paintComponent(g);
//
//            Graphics2D g2D = (Graphics2D) g;
//
//            if (scatterPlot != null) {
//
//                g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
//                        PreferenceManager.getInstance().getAntiAliasingHint());
//
//                Color color = g.getColor();
//                Font font = g.getFont();
//
//                Rectangle pointShape = scatterPlot.pointShape;
//
//                String categoryName = scatterPlot.dataModel.categoryName;
//                if (categoryName == null || categoryName.equals("")) {
//                    return;
//                }
//
//                g2D.setFont(headerFont);
//                g2D.drawString(categoryName, leftMargin, topMargin);
//
//                g2D.setFont(labelFont);
//                int y = topMargin + 20;
//                for (String sn : scatterPlot.dataModel.getSeriesNames()) {
//                    Color c = ScatterPlot.getColor(categoryName, sn);
//                    g2D.setColor(c);
//                    g2D.fillOval(leftMargin + 5, y, pointShape.width, pointShape.height);
//
//                    String displayString = sn.equals("") ? "Unknown" : sn;
//                    g2D.setColor(Color.black);
//                    g2D.drawString(displayString, leftMargin + 20, y + 5);
//
//                    y += 20;
//
//                }
//
//
//                g2D.setColor(color);
//                g2D.setFont(font);
//
//
//            }
//        }


    }

    static class LegendIcon implements Icon {

        Rectangle shape;
        Color color;

        LegendIcon(Rectangle shape, Color color) {
            this.shape = shape;
            this.color = color;
        }

        public void paintIcon(Component component, Graphics graphics, int x, int y) {
            Color c = graphics.getColor();
            graphics.setColor(color);
            graphics.fillOval(x, y, shape.width, shape.height);
            graphics.setColor(c);
        }

        public int getIconWidth() {
            return shape.width;
        }

        public int getIconHeight() {
            return shape.height;
        }
    }


//    class ToolBar extends JPanel {
//
//        ToolBar() {
//            init();
//        }
//
//        void init() {
//
//            setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
//            final JButton lassoButton = new JButton("Lasso");
//            lassoButton.addActionListener(new ActionListener() {
//                public void actionPerformed(ActionEvent actionEvent) {
//                    if (lassoInProgress == false) {
//                        lassoInProgress = true;
//                        lassoPath = new SelectionPath();
//                    } else {
//                        lassoInProgress = false;
//                        lassoPath = null;
//                    }
//                }
//            });
//            add(new JPanel(null));  // Spacer
//            add(lassoButton);
//
//        }
//    }

}
