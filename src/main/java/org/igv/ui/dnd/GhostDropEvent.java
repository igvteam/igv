package org.igv.ui.dnd;

import org.igv.track.Track;
import org.igv.ui.panel.TrackPanel;

import java.util.ArrayList;
import java.util.List;
import java.awt.Point;


public class GhostDropEvent {
    private final Point point;
    private final Point startPoint;
    List<Track> tracks;
    List<TrackPanel> sourcePanels;
    boolean tracksDropped = false;


    public GhostDropEvent(Point startPoint, Point point, List<Track> tracks) {
        this.startPoint = startPoint;
        this.point = point;
        this.tracks = tracks;
        sourcePanels = new ArrayList<>();
    }

    public void addSourcePanel(TrackPanel panel) {
        sourcePanels.add(panel);
    }

    public void setTracksDropped(boolean bool) {
        tracksDropped = bool;
    }

    public boolean isTracksDropped() {
        return tracksDropped;
    }

    /**
     * Called when the destination panel is found.  If the tracks are dropped outside a valid destination panel
     * this is never called
     */
    public void removeTracksFromSource() {

        for(TrackPanel panel : sourcePanels) {
            panel.removeTracks(tracks);
        }
        sourcePanels.clear();

    }

    public Point getStartLocation() {
        return startPoint;
    }

    public Point getDropLocation() {
        return point;
    }

    public java.util.List<Track> getTracks() {
        return tracks;
    }
}
