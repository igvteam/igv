package org.igv.ui.dnd;

import javax.swing.*;
import java.awt.*;


public abstract class AbstractGhostDropManager implements GhostDropListener {
    protected JComponent component;

    public AbstractGhostDropManager() {
        this(null);
    }

    public AbstractGhostDropManager(JComponent component) {
        this.component = component;
    }

    protected Point getTranslatedPoint(Point point) {
        Point p = (Point) point.clone();
        SwingUtilities.convertPointFromScreen(p, component);
        return p;
    }

    protected boolean isInTarget(Point point) {
        Rectangle bounds = component.getVisibleRect();
        return bounds.contains(point);
    }

    public void ghostDropped(GhostDropEvent e) {
    }
}
