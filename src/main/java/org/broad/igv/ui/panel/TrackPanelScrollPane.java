/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import com.jidesoft.swing.JideScrollPane;
import org.broad.igv.logging.*;
import org.broad.igv.ui.util.SnapshotUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

/**
 * @author jrobinso
 */
public class TrackPanelScrollPane extends JideScrollPane implements Paintable {

    private static Logger log = LogManager.getLogger(TrackPanelScrollPane.class);

    TrackPanel trackPanel;
    boolean isScrolling = false;

    public TrackPanelScrollPane() {
        setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(102, 102, 102)));
        setForeground(new java.awt.Color(153, 153, 153));
        setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        getVerticalScrollBar().setUnitIncrement(16);

        addMouseWheelListener(new MouseWheelListener() {
            public void mouseWheelMoved(MouseWheelEvent mouseWheelEvent) {
                trackPanel.getNamePanel().repaint();
            }
        });

        // A fix for name panel painting problems.   Not sure why this is neccessary, but it is.
        // The adustment listener forces a repaint after a scroll action is complete.
        final JScrollBar sb = this.getVerticalScrollBar();
        sb.addAdjustmentListener(new AdjustmentListener() {
            public void adjustmentValueChanged(AdjustmentEvent adjustmentEvent) {
                if (isScrolling) {
                    if (!adjustmentEvent.getValueIsAdjusting()) {
                        isScrolling = false;
                        trackPanel.getNamePanel().repaint();
                    }
                } else {
                    isScrolling = adjustmentEvent.getValueIsAdjusting();
                }

            }
        });
    }

    @Override
    public void setBackground(Color color) {
        super.setBackground(color);
        if (trackPanel != null)
            trackPanel.setBackground(color);
    }

    @Override
    public void setViewportView(Component trackSetView) {
        if (!(trackSetView instanceof TrackPanel)) {
            throw new IllegalArgumentException("Class TrackPanelScrollPane can only contain a TrackPanel");
        }
        super.setViewportView(trackSetView);
        this.trackPanel = (TrackPanel) trackSetView;
    }

    public TrackPanel getTrackPanel() {
        return trackPanel;
    }


    public String getTrackPanelName() {
        return trackPanel.getName();
    }

    public void minimizeHeight() {
        int prefHeight = trackPanel.getPreferredPanelHeight();
        if (prefHeight < trackPanel.getViewportHeight()) {
            this.setSize(getWidth(), prefHeight);
        }
    }

    public DataPanelContainer getDataPanel() {
        return trackPanel.getDataPanelContainer();
    }

    public boolean isEmpty() {
        return !trackPanel.hasTracks();
    }

    public TrackNamePanel getNamePanel() {
        return trackPanel.getNamePanel();
    }

    @Override
    public void paintOffscreen(Graphics2D g, Rectangle tspRect, boolean batch) {
        Graphics2D panelGraphics = (Graphics2D) g.create();
        int yOffset = getViewport().getViewPosition().y;
        panelGraphics.translate(0, -yOffset);
        Rectangle panelRect = getViewportBorderBounds();
        panelRect.y = yOffset;//new Rectangle(0, 0, trackPanel.getWidth(), getViewport().getHeight());
        panelRect.height = tspRect.height;    // Offscreen images can exceed viewport bounds
        trackPanel.paintOffscreen(panelGraphics, panelRect, batch);
        panelGraphics.dispose();
    }

    @Override
    public int getSnapshotHeight(boolean batch) {

        if(trackPanel.getTracks().size() == 0) {
            return 0;
        }

        if (batch) {
            int panelHeight;
            int maxPanelHeight = SnapshotUtilities.getMaxPanelHeight();
            final int scrollPaneHeight = getHeight();
            if (maxPanelHeight <= 0) {
                panelHeight = scrollPaneHeight;
            } else {
                int contentHeight = trackPanel.getPreferredPanelHeight();
                panelHeight = Math.min(maxPanelHeight, Math.max(scrollPaneHeight, contentHeight));
            }
            return panelHeight;
        } else {
            return getHeight();  // This is the height of the scroll pane, the visible height in the UI
        }
    }
}
