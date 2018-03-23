
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

import java.awt.*;


/**
 * Layout manager for chart panels, based loosely on BorderLayout
 */
public class ChartLayout implements LayoutManager2,
        java.io.Serializable {


    public static final String CHART = "Chart";
    public static final String XAXIS = "XAxis";
    public static final String YAXIS = "YAxis";
    public static final String TITLE = "Title";
    public static final String LEGEND = "Legend";

    int hgap = 2;
    int vgap = 2;


    Component title;
    Component yAxis;
    Component legend;
    Component xAxis;
    Component chart;

    public void addLayoutComponent(Component comp, Object constraints) {
        synchronized (comp.getTreeLock()) {
            if ((constraints == null) || (constraints instanceof String)) {
                addLayoutComponent((String) constraints, comp);
            } else {
                throw new IllegalArgumentException("cannot add to layout: constraint must be a string (or null)");
            }
        }
    }

    public void addLayoutComponent(String name, Component comp) {
        synchronized (comp.getTreeLock()) {
            if (CHART.equals(name)) {
                chart = comp;
            } else if (TITLE.equals(name)) {
                title = comp;
            } else if (XAXIS.equals(name)) {
                xAxis = comp;
            } else if (YAXIS.equals(name)) {
                yAxis = comp;
            } else if (LEGEND.equals(name)) {
                legend = comp;
            } else {
                throw new IllegalArgumentException("Unknown chart component: " + name);
            }
        }
    }

    /**
     * Removes the specified component from this border layout. This
     * method is called when a container calls its <code>remove</code> or
     * <code>removeAll</code> methods. Most applications do not call this
     * method directly.
     *
     * @param comp the component to be removed.
     * @see java.awt.Container#remove(java.awt.Component)
     * @see java.awt.Container#removeAll()
     */
    public void removeLayoutComponent(Component comp) {
        synchronized (comp.getTreeLock()) {
            if (comp == chart) {
                chart = null;
            } else if (comp == title) {
                title = null;
            } else if (comp == xAxis) {
                xAxis = null;
            } else if (comp == legend) {
                legend = null;
            } else if (comp == yAxis) {
                yAxis = null;
            }
        }
    }


    /**
     * Determines the minimum size of the <code>target</code> container
     * using this layout manager.
     * <p/>
     * This method is called when a container calls its
     * <code>getMinimumSize</code> method. Most applications do not call
     * this method directly.
     *
     * @param target the container in which to do the layout.
     * @return the minimum dimensions needed to lay out the subcomponents
     *         of the specified container.
     * @see java.awt.Container
     * @see java.awt.BorderLayout#preferredLayoutSize
     * @see java.awt.Container#getMinimumSize()
     */
    public Dimension minimumLayoutSize(Container target) {
        synchronized (target.getTreeLock()) {
            Dimension dim = new Dimension(0, 0);

            if (legend != null) {
                Dimension d = legend.getMinimumSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if (xAxis != null) {
                Dimension d = xAxis.getMinimumSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if (chart != null) {
                Dimension d = chart.getMinimumSize();
                dim.width += d.width;
                dim.height = Math.max(d.height, dim.height);
            }
            if (title != null) {
                Dimension d = title.getMinimumSize();
                //dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }
            if (yAxis != null) {
                Dimension d = yAxis.getMinimumSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }

            Insets insets = target.getInsets();
            dim.width += insets.left + insets.right;
            dim.height += insets.top + insets.bottom;

            return dim;
        }
    }

    /**
     * Determines the preferred size of the <code>target</code>
     * container using this layout manager, based on the components
     * in the container.
     * <p/>
     * Most applications do not call this method directly. This method
     * is called when a container calls its <code>getPreferredSize</code>
     * method.
     *
     * @param target the container in which to do the layout.
     * @return the preferred dimensions to lay out the subcomponents
     *         of the specified container.
     * @see java.awt.Container
     * @see java.awt.BorderLayout#minimumLayoutSize
     * @see java.awt.Container#getPreferredSize()
     */
    public Dimension preferredLayoutSize(Container target) {
        synchronized (target.getTreeLock()) {
            Dimension dim = new Dimension(0, 0);

            if (legend != null) {
                Dimension d = legend.getPreferredSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if (xAxis != null) {
                Dimension d = xAxis.getPreferredSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if (chart != null) {
                Dimension d = chart.getPreferredSize();
                dim.width += d.width;
                dim.height = Math.max(d.height, dim.height);
            }
            if (title != null) {
                Dimension d = title.getPreferredSize();
                //dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }
            if (yAxis != null) {
                Dimension d = yAxis.getPreferredSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }

            Insets insets = target.getInsets();
            dim.width += insets.left + insets.right;
            dim.height += insets.top + insets.bottom;

            return dim;
        }
    }

    /**
     * Returns the maximum dimensions for this layout given the components
     * in the specified target container.
     *
     * @param target the component which needs to be laid out
     * @see Container
     * @see #minimumLayoutSize
     * @see #preferredLayoutSize
     */
    public Dimension maximumLayoutSize(Container target) {
        return new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE);
    }

    /**
     * Returns the alignment along the x axis.  This specifies how
     * the component would like to be aligned relative to other
     * components.  The value should be a number between 0 and 1
     * where 0 represents alignment along the origin, 1 is aligned
     * the furthest away from the origin, 0.5 is centered, etc.
     */
    public float getLayoutAlignmentX(Container parent) {
        return 0.5f;
    }

    /**
     * Returns the alignment along the y axis.  This specifies how
     * the component would like to be aligned relative to other
     * components.  The value should be a number between 0 and 1
     * where 0 represents alignment along the origin, 1 is aligned
     * the furthest away from the origin, 0.5 is centered, etc.
     */
    public float getLayoutAlignmentY(Container parent) {
        return 0.5f;
    }

    /**
     * Invalidates the layout, indicating that if the layout manager
     * has cached information it should be discarded.
     */
    public void invalidateLayout(Container target) {
    }

    /**
     * Lays out the container argument. The preferred width and height of the components are used selectively
     * as follows
     *
     * Heights:  title, xAxis
     * Widths:   yAxis, legend
     *
     * The chart gets whatever is left over.
     *
     *
     * @param target the container in which to do the layout.
     * @see java.awt.Container
     * @see java.awt.Container#doLayout()
     */
    public void layoutContainer(Container target) {
        synchronized (target.getTreeLock()) {
            Insets insets = target.getInsets();
            int top = insets.top;
            int bottom = target.getHeight() - insets.bottom;
            int left = insets.left;
            int right = target.getWidth() - insets.right;

            int titleHeight = title == null ? 0 : title.getPreferredSize().height;
            int xAxisHeight = xAxis == null ? 0 : xAxis.getPreferredSize().height;

            int yAxisWidth = yAxis == null ? 0 : yAxis.getPreferredSize().width;
            int legendWidth = legend == null ? 0 : legend.getPreferredSize().width;

            int x1 = left + yAxisWidth +  (yAxisWidth == 0 ? 0 : hgap);
            int x2 = right - legendWidth - (legendWidth == 0 ? 0 : hgap);

            int y1 = top + titleHeight + (titleHeight == 0 ? 0 : vgap);
            int y2 = bottom - xAxisHeight - (xAxisHeight == 0 ? 0 : vgap);


            if (title != null) {
                title.setBounds(left, top, right - left, titleHeight);
            }
            if (xAxis != null) {
              //   xAxis.setSize(new Dimension(x2-x1, xAxisHeight));
              //  xAxis.setSize(x2-x1, xAxisHeight);
                xAxis.setBounds(x1, y2 + vgap, x2 - x1, xAxisHeight);
            }
            if (legend != null) {
                legend.setBounds(right - legendWidth, y1, legendWidth, y2 - y1);
            }
            if (yAxis != null) {
                yAxis.setBounds(left, y1, yAxisWidth, y2 - y1);
            }
            if (chart != null) {
              //  chart.setSize(new Dimension(x2-x1, y2-y1));
              //  chart.setSize(x2-x1, y2-y1);
                chart.setBounds(x1, y1, x2 - x1, y2 - y1);
            }
        }
    }

}
