package org.igv.ui.panel;

import org.igv.event.IGVEventBus;
import org.igv.event.TrackSelectionEvent;
import org.igv.track.Track;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.List;

/**
 * A panel that contains a checkbox for selecting a track.
 * This panel is displayed to the left of the drag handle in TrackPanel
 * and is invisible by default. It can be activated from a menu item.
 * The checkbox state is the source of truth for track selection.
 */
public class TrackSelectionPanel extends JPanel {

    public static final int SELECTION_PANEL_WIDTH = 24;

    /** Panel of the last ordinary (non-Shift) selection click; one end of a Shift-click range. */
    private static TrackSelectionPanel anchor;

    private final TrackPanel trackPanel;
    private final JCheckBox checkBox;

    public TrackSelectionPanel(TrackPanel trackPanel) {
        this.trackPanel = trackPanel;
        setBackground(UIConstants.getTrackPanelBackground());
        setPreferredSize(new Dimension(SELECTION_PANEL_WIDTH, 0));
        setMinimumSize(new Dimension(SELECTION_PANEL_WIDTH, 0));
        setLayout(new GridBagLayout());

        checkBox = new JCheckBox();
        checkBox.setBackground(getBackground());
        checkBox.setOpaque(true);
        checkBox.addItemListener(e -> IGVEventBus.getInstance().post(new TrackSelectionEvent()));
        // Action events fire only for user clicks, not setSelected(), so programmatic changes never move the anchor
        checkBox.addActionListener(e -> {
            boolean shift = (e.getModifiers() & ActionEvent.SHIFT_MASK) != 0;
            if (!(shift && selectRangeFromAnchor())) {
                anchor = this;
            }
        });

        add(checkBox);

        // Initially invisible
        setVisible(false);
    }

    /**
     * Toggle this track's selection in response to a user click.  With {@code extendRange} (Shift-click),
     * instead select every visible track between the last ordinary click and this one.
     */
    public void toggleTrackSelection(boolean extendRange) {
        if (extendRange && selectRangeFromAnchor()) {
            return;
        }
        setTrackSelected(!isTrackSelected());
        anchor = this;
    }

    /**
     * Select the inclusive range of tracks, in visual order, from the anchor to this track.
     *
     * @return false if there is no usable anchor, in which case nothing is changed
     */
    private boolean selectRangeFromAnchor() {
        if (anchor == null || anchor == this) {
            return false;
        }
        List<TrackPanel> trackPanels = IGV.getInstance().getMainPanel().getTrackPanels();
        int anchorIndex = trackPanels.indexOf(anchor.trackPanel);
        int index = trackPanels.indexOf(trackPanel);
        if (anchorIndex < 0 || index < 0) {
            return false;
        }
        for (int i = Math.min(anchorIndex, index); i <= Math.max(anchorIndex, index); i++) {
            TrackPanel tp = trackPanels.get(i);
            TrackPanelScrollPane sp = tp.getScrollPane();
            if (sp != null && sp.getSelectionPanel() != null && tp.getTrack().isVisible()) {
                sp.getSelectionPanel().setTrackSelected(true);
            }
        }
        return true;
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
