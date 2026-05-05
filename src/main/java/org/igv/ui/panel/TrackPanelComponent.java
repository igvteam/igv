package org.igv.ui.panel;


import org.igv.Globals;
import org.igv.track.Track;
import org.igv.track.TrackClickEvent;
import org.igv.track.TrackMenuUtils;
import org.igv.ui.IGV;

import javax.swing.*;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;


/**
 * @author eflakes
 */
abstract public class TrackPanelComponent extends JPanel {

    protected final boolean darkMode;
    List<MouseableRegion> mouseRegions;
    private TrackPanel trackPanel;

    public TrackPanelComponent(TrackPanel trackPanel) {
        this.trackPanel = trackPanel;
        setFocusable(true);
        mouseRegions = new ArrayList<>();
        this.darkMode = Globals.isDarkMode();
    }

    public TrackPanel getTrackPanel() {
        return trackPanel;
    }

    public Track getTrack() {
        return getTrackPanel().getTrack();
    }

    public String getTrackSetID() {
        return getTrackPanel().getName();
    }

    protected void removeMousableRegions() {
        mouseRegions.clear();
    }

    protected List<MouseableRegion> getMouseRegions() {
        return mouseRegions;
    }

    public List<Track> getAllTracks() {
        TrackPanel dataTrackView = (TrackPanel) getParent();
        return dataTrackView.getTracks();
    }

    protected void openPopupMenu(TrackClickEvent te) {

        MouseEvent e = te.getMouseEvent();

        Track track = getTrack();

        IGVPopupMenu menu = TrackMenuUtils.getPopupMenu(track, track.getName(), te);

        if (menu != null) {
            menu.show(e.getComponent(), e.getX(), e.getY());
        }

    }


    public void saveImage(String extension) {
        IGV.getInstance().saveImage(getTrackPanel().getScrollPane(), "igv_panel", extension);
    }


}
