package org.igv.ui.panel;

import org.igv.event.IGVEventBus;
import org.igv.event.TrackSelectionEvent;
import org.igv.track.Track;
import org.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;

/**
 * A panel that contains a checkbox for selecting a track.
 * This panel is displayed to the left of the drag handle in TrackPanel
 * and is invisible by default. It can be activated from a menu item.
 * The checkbox state is the source of truth for track selection.
 */
public class TrackSelectionPanel extends JPanel {

    public static final int SELECTION_PANEL_WIDTH = 24;

    // Static flag to control visibility of all selection panels
    private static boolean selectionModeActive = false;

    private final TrackPanel trackPanel;
    private final JCheckBox checkBox;

    public TrackSelectionPanel(TrackPanel trackPanel) {
        this.trackPanel = trackPanel;
        setBackground(Color.WHITE);
        setPreferredSize(new Dimension(SELECTION_PANEL_WIDTH, 0));
        setMinimumSize(new Dimension(SELECTION_PANEL_WIDTH, 0));
        setLayout(new GridBagLayout());

        checkBox = new JCheckBox();
        checkBox.setBackground(Color.WHITE);
        checkBox.setOpaque(true);
        checkBox.addItemListener(e -> IGVEventBus.getInstance().post(new TrackSelectionEvent()));

        add(checkBox);

        // Initially invisible
        setVisible(false);
    }

    /**
     * Check if the track is selected (checkbox is checked)
     */
    public boolean isTrackSelected() {
        return checkBox.isSelected();
    }

    /**
     * Set the selection state of this track
     */
    public void setTrackSelected(boolean selected) {
        checkBox.setSelected(selected);
    }

    /**
     * Get the track associated with this selection panel
     */
    public Track getTrack() {
        return trackPanel.getTrack();
    }

    /**
     * Check if selection mode is active globally
     */
    public static boolean isSelectionModeActive() {
        return selectionModeActive;
    }

    /**
     * Set the global selection mode
     * @param active true to show selection checkboxes, false to hide them
     */
    public static void setSelectionModeActive(boolean active) {
        selectionModeActive = active;
    }

    /**
     * Toggle the selection mode
     * @return the new state of selection mode
     */
    public static boolean toggleSelectionMode() {
        selectionModeActive = !selectionModeActive;
        return selectionModeActive;
    }

    /**
     * Get the width of the selection panel (0 if not visible)
     */
    public int getEffectiveWidth() {
        return isVisible() ? SELECTION_PANEL_WIDTH : 0;
    }

    @Override
    public void setBackground(Color bg) {
        super.setBackground(bg);
        if (checkBox != null) {
            checkBox.setBackground(bg);
        }
    }

    public JCheckBox getCheckBox() {
        return checkBox;
    }

    public TrackPanel getTrackPanel() {
        return trackPanel;
    }
}
