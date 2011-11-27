package org.broad.igv.hic;

import java.awt.*;


/**
 * Layout manager for HiC panels, based loosely on BorderLayout
 */
public class HiCLayout implements LayoutManager2,
        java.io.Serializable {
    int hgap;
    int vgap;

    Component north;
    Component west;
    Component east;
    Component south;
    Component center;


    /**
     * The north layout constraint (top of container).
     */
    public static final String NORTH = "North";

    /**
     * The south layout constraint (bottom of container).
     */
    public static final String SOUTH = "South";

    /**
     * The east layout constraint (right side of container).
     */
    public static final String EAST = "East";

    /**
     * The west layout constraint (left side of container).
     */
    public static final String WEST = "West";

    /**
     * The center layout constraint (middle of container).
     */
    public static final String CENTER = "Center";


    /**
     * Constructs a new border layout with
     * no gaps between components.
     */
    public HiCLayout() {
        this(0, 0);
    }

    /**
     * Constructs a border layout with the specified gaps
     * between components.
     * The horizontal gap is specified by <code>hgap</code>
     * and the vertical gap is specified by <code>vgap</code>.
     *
     * @param hgap the horizontal gap.
     * @param vgap the vertical gap.
     */
    public HiCLayout(int hgap, int vgap) {
        this.hgap = hgap;
        this.vgap = vgap;
    }

    /**
     * Returns the horizontal gap between components.
     *
     * @since JDK1.1
     */
    public int getHgap() {
        return hgap;
    }

    /**
     * Sets the horizontal gap between components.
     *
     * @param hgap the horizontal gap between components
     * @since JDK1.1
     */
    public void setHgap(int hgap) {
        this.hgap = hgap;
    }

    /**
     * Returns the vertical gap between components.
     *
     * @since JDK1.1
     */
    public int getVgap() {
        return vgap;
    }

    /**
     * Sets the vertical gap between components.
     *
     * @param vgap the vertical gap between components
     * @since JDK1.1
     */
    public void setVgap(int vgap) {
        this.vgap = vgap;
    }

    /**
     * Adds the specified component to the layout, using the specified
     * constraint object.  For border layouts, the constraint must be
     * one of the following constants:  <code>NORTH</code>,
     * <code>SOUTH</code>, <code>EAST</code>,
     * <code>WEST</code>, or <code>CENTER</code>.
     * <p/>
     * Most applications do not call this method directly. This method
     * is called when a component is added to a container using the
     * <code>Container.add</code> method with the same argument types.
     *
     * @param comp        the component to be added.
     * @param constraints an object that specifies how and where
     *                    the component is added to the layout.
     * @throws IllegalArgumentException if the constraint object is not
     *                                  a string, or if it not one of the five specified
     *                                  constants.
     * @see java.awt.Container#add(java.awt.Component, java.lang.Object)
     * @since JDK1.1
     */
    public void addLayoutComponent(Component comp, Object constraints) {
        synchronized (comp.getTreeLock()) {
            if ((constraints == null) || (constraints instanceof String)) {
                addLayoutComponent((String) constraints, comp);
            } else {
                throw new IllegalArgumentException("cannot add to layout: constraint must be a string (or null)");
            }
        }
    }

    /**
     * @deprecated replaced by <code>addLayoutComponent(Component, Object)</code>.
     */
    @Deprecated
    public void addLayoutComponent(String name, Component comp) {
        synchronized (comp.getTreeLock()) {
            /* Special case:  treat null the same as "Center". */
            if (name == null) {
                name = "Center";
            }

            /* Assign the component to one of the known regions of the layout.
            */
            if ("Center".equals(name)) {
                center = comp;
            } else if ("North".equals(name)) {
                north = comp;
            } else if ("South".equals(name)) {
                south = comp;
            } else if ("East".equals(name)) {
                east = comp;
            } else if ("West".equals(name)) {
                west = comp;
            } else {
                throw new IllegalArgumentException("cannot add to layout: unknown constraint: " + name);
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
            if (comp == center) {
                center = null;
            } else if (comp == north) {
                north = null;
            } else if (comp == south) {
                south = null;
            } else if (comp == east) {
                east = null;
            } else if (comp == west) {
                west = null;
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

            boolean ltr = target.getComponentOrientation().isLeftToRight();
            Component c = null;

            if ((c = getChild(EAST, ltr)) != null) {
                Dimension d = c.getMinimumSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((c = getChild(WEST, ltr)) != null) {
                Dimension d = c.getMinimumSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((c = getChild(CENTER, ltr)) != null) {
                Dimension d = c.getMinimumSize();
                dim.width += d.width;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((c = getChild(NORTH, ltr)) != null) {
                Dimension d = c.getMinimumSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }
            if ((c = getChild(SOUTH, ltr)) != null) {
                Dimension d = c.getMinimumSize();
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

            boolean ltr = target.getComponentOrientation().isLeftToRight();
            Component c = null;

            if ((c = getChild(EAST, ltr)) != null) {
                Dimension d = c.getPreferredSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((c = getChild(WEST, ltr)) != null) {
                Dimension d = c.getPreferredSize();
                dim.width += d.width + hgap;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((c = getChild(CENTER, ltr)) != null) {
                Dimension d = c.getPreferredSize();
                dim.width += d.width;
                dim.height = Math.max(d.height, dim.height);
            }
            if ((c = getChild(NORTH, ltr)) != null) {
                Dimension d = c.getPreferredSize();
                dim.width = Math.max(d.width, dim.width);
                dim.height += d.height + vgap;
            }
            if ((c = getChild(SOUTH, ltr)) != null) {
                Dimension d = c.getPreferredSize();
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
     * Lays out the container argument using this border layout.
     * <p/>
     * This method actually reshapes the components in the specified
     * container in order to satisfy the constraints of this
     * <code>BorderLayout</code> object. The <code>NORTH</code>
     * and <code>SOUTH</code> components, if any, are placed at
     * the top and bottom of the container, respectively. The
     * <code>WEST</code> and <code>EAST</code> components are
     * then placed on the left and right, respectively. Finally,
     * the <code>CENTER</code> object is placed in any remaining
     * space in the middle.
     * <p/>
     * Most applications do not call this method directly. This method
     * is called when a container calls its <code>doLayout</code> method.
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

            int dh = bottom - top;
            int dw = right - left;


            // reserve space for chr1 & chr2 legneds (north and west)
            int dx = west == null ? 0 : west.getPreferredSize().width + hgap;
            int dy = north == null ? 0 : north.getPreferredSize().height + vgap;
            if (south != null) dy += south.getPreferredSize().height + vgap;

            // Center must be square
            int centerDim = Math.min(dw - dx, dh - dy);

            if (north != null) {
                Dimension d = north.getPreferredSize();
                north.setBounds(dx, top, centerDim, d.height);
                top += d.height + vgap;
            }
            if (west != null) {
                Dimension d = west.getPreferredSize();
                west.setBounds(left, top, d.width, centerDim);
            }
            if (center != null) {
                center.setBounds(left + west.getWidth() + hgap, top, centerDim, centerDim);
            }
            if (east != null) {
                int xEast = left + west.getWidth() + hgap + centerDim + hgap;
                int wEast = Math.max(0, dw - xEast);
                east.setBounds(xEast, top, wEast, centerDim);
            }
            if (south != null) {
                top += centerDim + vgap;
                int hSouth = Math.max(0, dh - top);
                south.setBounds(left, top, dw, hSouth);
            }

        }
    }

    /**
     * Get the component that corresponds to the given constraint location
     *
     * @param key The desired absolute position,
     *            either NORTH, SOUTH, EAST, or WEST.
     * @param ltr Is the component line direction left-to-right?
     */
    private Component getChild(String key, boolean ltr) {
        Component result = null;

        if (key == NORTH) {
            result = north;
        } else if (key == SOUTH) {
            result = south;
        } else if (key == WEST) {

            result = west;

        } else if (key == EAST) {
            result = east;

        } else if (key == CENTER) {
            result = center;
        }
        if (result != null && !result.isVisible()) {
            result = null;
        }
        return result;
    }

    /**
     * Returns a string representation of the state of this border layout.
     *
     * @return a string representation of this border layout.
     */
    public String toString() {
        return getClass().getName() + "[hgap=" + hgap + ",vgap=" + vgap + "]";
    }
}
