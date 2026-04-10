package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;

import javax.swing.*;
import java.awt.*;

/**
 * A scrollable panel container that uses BoxLayout and doesn't stretch vertically
 * when placed in a JScrollPane viewport.
 */
public class ScrollableTrackContainer extends JPanel implements Scrollable {

    private MainPanel mainPanel;

    public ScrollableTrackContainer(MainPanel mainPanel) {
        this.mainPanel = mainPanel;
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        if (Globals.isDarkMode() && !PreferencesManager.getPreferences().hasExplicitValue(Constants.BACKGROUND_COLOR)) {
            setBackground(UIManager.getColor("Panel.background"));
        } else {
            setBackground(PreferencesManager.getPreferences().getAsColor(Constants.BACKGROUND_COLOR));
        }
    }

    @Override
    protected void paintChildren(Graphics g) {
        super.paintChildren(g);
        drawPanelDividers(g);
    }

    /**
     * Draw vertical divider lines at the boundaries between name, attribute, and data panels.
     * Uses computed positions from MainPanel so lines update immediately after layout changes
     * (e.g. hiding the attribute panel or changing name panel width).
     */
    private void drawPanelDividers(Graphics g) {
        int leftOffset = mainPanel.getLeftOffset();
        int nameRight = leftOffset + mainPanel.getNamePanelX() + mainPanel.getNamePanelWidth();
        int dataLeft = leftOffset + mainPanel.getDataPanelX();

        int h = getHeight();
        Color dividerColor = Globals.isDarkMode() ? Color.GRAY : Color.LIGHT_GRAY;
        g.setColor(dividerColor);

        g.drawLine(nameRight, 0, nameRight, h);
        if (dataLeft != nameRight) {
            g.drawLine(dataLeft, 0, dataLeft, h);
        }
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

