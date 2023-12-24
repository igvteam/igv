package org.broad.igv.ui.panel;

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import javax.swing.*;
import java.awt.*;
import java.util.logging.Level;

public class MainPanelLayout implements LayoutManager2 {

    private static Logger log = LogManager.getLogger(MainPanelLayout.class);
    Container target;

    int minSize = -1;
    int prefSize = -1;
    int maxSize = -1;

    int width = -1;

    public MainPanelLayout(Container target) {
        this.target = target;
    }

    @Override
    /**
     * Adds the specified component to the layout, using the specified
     * constraint object.  Called from the owning Container
     *
     * @param comp the component to be added
     * @param constraints  where/how the component is added to the layout.
     */
    public void addLayoutComponent(Component comp, Object constraints) {
        invalidateLayout(comp.getParent());
    }

    /**
     * Calculates the maximum size dimensions for the specified container,
     * given the components it contains.
     *
     * @param target the target container
     * @return the maximum size of the container
     * @see java.awt.Component#getMaximumSize
     * @see LayoutManager
     */
    @Override
    public Dimension maximumLayoutSize(Container target) {
        checkSizes();
        return new Dimension(width, maxSize);
    }

    /**
     * Returns the alignment along the x axis.  This specifies how the component would like to be aligned relative to
     * other components.  The value should be a number between 0 and 1 where 0 represents alignment along the origin, 1
     * is aligned the furthest away from the origin, 0.5 is centered, etc.
     */
    @Override
    public float getLayoutAlignmentX(Container target) {
        return 0;
    }

    /**
     * Returns the alignment along the y axis.  This specifies how the component would like to be aligned relative to
     * other components.  The value should be a number between 0 and 1 where 0 represents alignment along the origin, 1
     * is aligned the furthest away from the origin, 0.5 is centered, etc.
     */
    @Override
    public float getLayoutAlignmentY(Container target) {
        return 0;
    }

    /**
     * Invalidates the layout, indicating that if the layout manager
     * has cached information it should be discarded.
     *
     * @param target the target container
     */
    @Override
    public void invalidateLayout(Container target) {
        minSize = -1;
    }

    @Override
    public void addLayoutComponent(String name, Component comp) {
        invalidateLayout(comp.getParent());
    }

    @Override
    public void removeLayoutComponent(Component comp) {
        invalidateLayout(comp.getParent());
    }

    @Override
    public Dimension preferredLayoutSize(Container parent) {
        checkSizes();
        return new Dimension(width, prefSize);
    }

    @Override
    public Dimension minimumLayoutSize(Container parent) {
        checkSizes();
        return new Dimension(width, minSize);
    }

    @Override
    public void layoutContainer(Container parent) {
        int n = target.getComponentCount();

        int w = target.getWidth();
        int y = 0;
        for (int i = 0; i < n; i++) {
            Component c = target.getComponent(i);
            if (c.isVisible()) {

                Dimension d;
                if (c instanceof TrackPanelScrollPane) {
                    TrackPanelScrollPane sp = (TrackPanelScrollPane) c;
                    Insets insets = sp.getInsets();
                    int h = ((TrackPanelScrollPane) c).getTrackPanel().getTotalTrackHeight() + insets.bottom + insets.top;
                    d = new Dimension(w, h);
                } else {
                    // This should never happen
                    d = c.getPreferredSize();
                }

                c.setSize(d);
                c.setLocation(0, y);
                y += d.height;
            } else {
                c.setSize(new Dimension(0, 0));
            }
        }

        Dimension targetPrefSize = target.getPreferredSize();
        targetPrefSize.height = y;
        target.setPreferredSize(targetPrefSize);
    }

    void checkSizes() {
        if (minSize < 0) {
            minSize = 0;
            prefSize = 0;
            maxSize = 0;
            int n = target.getComponentCount();
            for (int i = 0; i < n; i++) {
                Component c = target.getComponent(i);
                if (c.isVisible()) {
                    minSize += c.getMinimumSize().height;
                    prefSize += c.getPreferredSize().height;
                    maxSize += c.getMaximumSize().height;
                }
            }
            minSize = Math.min(32767, minSize);
            prefSize = Math.min(32767, prefSize);
            maxSize = Math.min(32767, maxSize);
            width = target.getWidth();
        }
    }
}
