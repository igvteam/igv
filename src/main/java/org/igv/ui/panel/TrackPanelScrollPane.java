package org.igv.ui.panel;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ui.util.SnapshotUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;

/**
 * @author jrobinso
 */
public class TrackPanelScrollPane extends JScrollPane implements Paintable {

    private static Logger log = LogManager.getLogger(TrackPanelScrollPane.class);

    TrackPanel trackPanel;

    public TrackPanelScrollPane() {
        setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(102, 102, 102)));
        setForeground(new java.awt.Color(153, 153, 153));
        setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
        getVerticalScrollBar().setUnitIncrement(16);

        // Repaint name panels on end of scroll to center the text in the visible window.
        getVerticalScrollBar().addAdjustmentListener(new AdjustmentListener() {
            @Override
            public void adjustmentValueChanged(AdjustmentEvent e) {
                if (!e.getValueIsAdjusting()) {
                    trackPanel.getNamePanel().repaint();
                }
            }
        });
    }

    /**
     * Set an explicit height for this scroll pane. This overrides the content-based sizing
     * and locks the scroll pane to this height until changed by another call to this method.
     *
     */
    public void validateTrackHeight() {

        // Update scrollbar visibility based on new height
        updateScrollbarPolicy();

        // Invalidate ourselves first
        invalidate();

        // Revalidate the parent container to trigger BoxLayout recalculation
        // This will resize us to the new explicit height
        Container parent = getParent();
        if (parent != null) {
            parent.invalidate();
            // Get the root container to force a complete layout pass
            Component root = SwingUtilities.getRoot(this);
            if (root instanceof Container) {
                ((Container) root).validate();
            }
        }

        repaint();
    }

    /**
     * Updates the vertical scrollbar policy based on whether the content height
     * exceeds the scroll pane's preferred height. Call this when the track's
     * content height may have changed.
     */
    public void updateScrollbarPolicy() {
        if (trackPanel != null && trackPanel.getTrack() != null) {
            int contentHeight = trackPanel.getTrack().getContentHeight();
            int scrollPaneHeight = getPreferredSize().height;
            boolean needsScrolling = contentHeight > scrollPaneHeight;

            int newPolicy = needsScrolling
                    ? ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS
                    : ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER;

            if (getVerticalScrollBarPolicy() != newPolicy) {
                setVerticalScrollBarPolicy(newPolicy);
            }
        }
    }

    @Override
    public Dimension getPreferredSize() {
        int height = trackPanel.getTrack().getHeight();
        return new Dimension(super.getPreferredSize().width, height);
    }

    @Override
    public Dimension getMaximumSize() {
        // Limit maximum height to preferred height so BoxLayout doesn't stretch
        Dimension pref = getPreferredSize();
        return new Dimension(Integer.MAX_VALUE, pref.height);
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
        updateScrollbarPolicy();
    }

    public TrackPanel getTrackPanel() {
        return trackPanel;
    }

    public String getTrackPanelName() {
        return trackPanel.getName();
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

        if (trackPanel.getTracks().size() == 0) {
            return 0;
        }

        if (batch) {
            int panelHeight;
            int maxPanelHeight = SnapshotUtilities.getMaxPanelHeight();
            final int scrollPaneHeight = getHeight();
            if (maxPanelHeight <= 0) {
                panelHeight = scrollPaneHeight;
            } else {
                int contentHeight = trackPanel.getContentHeight();
                panelHeight = Math.min(maxPanelHeight, Math.max(scrollPaneHeight, contentHeight));
            }
            return panelHeight;
        } else {
            return getHeight();  // This is the height of the scroll pane, the visible height in the UI
        }
    }
}
