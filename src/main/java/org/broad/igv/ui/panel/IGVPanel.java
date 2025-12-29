/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.broad.igv.logging.*;

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

        TrackPanelScrollPane scollpane = null;
        Container parent = getParent();

        if (parent instanceof JViewport) {
            scollpane = (TrackPanelScrollPane) ((JViewport) parent).getParent();
        }
        return scollpane;
    }

    @Override
    public void doLayout() {

        synchronized (getTreeLock()) {

            int h = getHeight(); //getPreferredSize().height;
            Component[] children = getComponents();

            int nw = mainPanel.getNamePanelWidth();

            Component namePanel = children[0];
            Component attributePanel = children[1];
            Component dataPanel = children[2];
            namePanel.setBounds(mainPanel.getNamePanelX(), 0, nw, h);

            attributePanel.setBounds(mainPanel.getAttributePanelX(), 0, mainPanel.getAttributePanelWidth(), h);  // Attributes
            dataPanel.setBounds(mainPanel.getDataPanelX(), 0, mainPanel.getDataPanelWidth(), h);

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
