package org.igv.ui.panel;

import org.igv.event.IGVEvent;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
import org.igv.event.TrackSelectionEvent;
import org.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * Application-header counterpart to {@link TrackSelectionPanel}: a select-all
 * checkbox positioned in the leftmost {@link TrackSelectionPanel#SELECTION_PANEL_WIDTH}
 * region so it aligns with the per-track selection checkboxes. Toggling it
 * sets every track's selection state to match. The checkbox auto-updates to
 * reflect whether all tracks are currently selected.
 */
public class HeaderSelectAllPanel extends JPanel implements IGVEventObserver {

    private final JCheckBox checkBox;
    private boolean syncingFromTracks;

    public HeaderSelectAllPanel() {
        setLayout(null);
        setOpaque(false);

        checkBox = new JCheckBox();
        checkBox.setOpaque(false);
        checkBox.setToolTipText("Select / unselect all tracks");
        checkBox.addItemListener(e -> {
            if (syncingFromTracks) return;
            boolean selected = checkBox.isSelected();
            for (TrackPanel tp : IGV.getInstance().getMainPanel().getTrackPanels()) {
                TrackPanelScrollPane sp = tp.getScrollPane();
                if (sp != null && sp.getSelectionPanel() != null) {
                    sp.getSelectionPanel().setTrackSelected(selected);
                }
            }
        });
        add(checkBox);

        IGVEventBus.getInstance().subscribe(TrackSelectionEvent.class, this);
    }

    private static final int BOTTOM_MARGIN = 5;

    @Override
    public void doLayout() {
        int h = getHeight();
        int selW = TrackSelectionPanel.SELECTION_PANEL_WIDTH;
        Dimension cbSize = checkBox.getPreferredSize();
        int x = Math.max(0, (selW - cbSize.width) / 2);
        int y = Math.max(0, h - cbSize.height - BOTTOM_MARGIN);
        checkBox.setBounds(x, y, cbSize.width, cbSize.height);
    }

    @Override
    public void receiveEvent(IGVEvent event) {
        if (event instanceof TrackSelectionEvent) {
            syncStateFromTracks();
        }
    }

    private void syncStateFromTracks() {
        List<TrackPanel> trackPanels = IGV.getInstance().getMainPanel().getTrackPanels();
        boolean allSelected = !trackPanels.isEmpty();
        for (TrackPanel tp : trackPanels) {
            TrackPanelScrollPane sp = tp.getScrollPane();
            if (sp == null || sp.getSelectionPanel() == null || !sp.getSelectionPanel().isTrackSelected()) {
                allSelected = false;
                break;
            }
        }
        if (allSelected != checkBox.isSelected()) {
            syncingFromTracks = true;
            try {
                checkBox.setSelected(allSelected);
            } finally {
                syncingFromTracks = false;
            }
        }
    }

    /**
     * Show or hide the select-all checkbox. Called when the
     * {@code SHOW_SELECTION_PANEL} preference is toggled.
     */
    public void setCheckBoxVisible(boolean visible) {
        if (!visible && checkBox.isSelected()) {
            syncingFromTracks = true;
            try {
                checkBox.setSelected(false);
            } finally {
                syncingFromTracks = false;
            }
        }
        checkBox.setVisible(visible);
    }
}
