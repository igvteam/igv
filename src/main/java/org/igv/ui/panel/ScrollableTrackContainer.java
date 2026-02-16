package org.igv.ui.panel;

import javax.swing.*;
import java.awt.*;

/**
 * A scrollable panel container that uses BoxLayout and doesn't stretch vertically
 * when placed in a JScrollPane viewport.
 */
public class ScrollableTrackContainer extends JPanel implements Scrollable {

    public ScrollableTrackContainer() {
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
    }

    // Scrollable implementation

    @Override
    public Dimension getPreferredSize() {
        // Calculate preferred size based on children
        int height = 0;
        int width = 0;
        for (Component c : getComponents()) {
            Dimension pref = c.getPreferredSize();
            height += pref.height;
            width = Math.max(width, pref.width);
        }
        return new Dimension(width, height);
    }

    @Override
    public Dimension getPreferredScrollableViewportSize() {
        return getPreferredSize();
    }

    @Override
    public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
        return 16;
    }

    @Override
    public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
        return orientation == SwingConstants.VERTICAL ? visibleRect.height : visibleRect.width;
    }

    @Override
    public boolean getScrollableTracksViewportWidth() {
        return true; // Stretch horizontally to fit viewport
    }

    @Override
    public boolean getScrollableTracksViewportHeight() {
        return false; // Do NOT stretch vertically - use preferred height
    }
}

