package org.igv.ui.panel;

import org.igv.Globals;
import org.igv.logging.*;

import javax.swing.*;
import java.awt.*;

/**
 * A panel class that lays out its components in a very specific way so that
 * children (i.e. name, attribute, and data panels) from all instances align vertically
 */
public class IGVPanel extends JPanel implements Paintable {

    Logger log = LogManager.getLogger(IGVPanel.class);

    // Backpointer to parent
    MainPanel mainPanel;

    public IGVPanel(MainPanel mainPanel) {
        setLayout(null); //new TrackPanelLayout());
        this.mainPanel = mainPanel;
    }

    public int getViewportHeight() {

        Container parent = getParent();
        return parent == null ? 0 : parent.getHeight();
    }

    public TrackPanelScrollPane getScrollPane() {
        Container parent = getParent();
        while (parent != null) {
            if (parent instanceof TrackPanelScrollPane) {
                return (TrackPanelScrollPane) parent;
            }
            parent = parent.getParent();
        }
        return null;
    }

    @Override
    protected void paintChildren(Graphics g) {
        super.paintChildren(g);
        // Only draw dividers on the header IGVPanel, not on TrackPanel subclasses.
        // TrackPanel positions children without leftOffset (drag handle is outside),
        // and track area dividers are drawn by ScrollableTrackContainer instead.
        if (!(this instanceof TrackPanel)) {
            drawPanelDividers(g);
        }
    }

    /**
     * Draw vertical divider lines at the boundaries between name, attribute, and data panels.
     * Uses the same coordinate computation as ScrollableTrackContainer for exact alignment.
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


    @Override
    public void doLayout() {

        synchronized (getTreeLock()) {

            int h = getHeight(); //getPreferredSize().height;
            Component[] children = getComponents();

            int leftOffset = mainPanel.getLeftOffset();
            int nw = mainPanel.getNamePanelWidth();

            Component dragHandleSpacer = children[0];
            Component namePanel = children[1];
            Component attributePanel = children[2];
            Component dataPanel = children[3];

            dragHandleSpacer.setBounds(0, 0, leftOffset, h);
            namePanel.setBounds(mainPanel.getNamePanelX() + leftOffset, 0, nw, h);
            attributePanel.setBounds(mainPanel.getAttributePanelX() + leftOffset, 0, mainPanel.getAttributePanelWidth(), h);
            dataPanel.setBounds(mainPanel.getDataPanelX() + leftOffset, 0, mainPanel.getDataPanelWidth() - leftOffset, h);

            dataPanel.doLayout();
        }
    }

    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        g.setColor(Color.black);

        Component[] children = getComponents();

        for(Component component : children) {
            if (component instanceof Paintable) {
                if (component.getWidth() > 0) {
                    Rectangle clipRect = new Rectangle(0, rect.y, component.getWidth(), rect.height);
                    Graphics2D graphics = (Graphics2D) g.create();
                    graphics.translate(component.getX(), component.getY());
                    //  graphics.setClip(clipRect);
                    ((Paintable) component).paintOffscreen(graphics, clipRect, batch);
                    graphics.dispose();
                }
            }
        }
    }

    @Override
    public int getSnapshotHeight(boolean batch) {
        return getHeight();
    }
}
