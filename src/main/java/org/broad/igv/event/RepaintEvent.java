package org.broad.igv.event;

import org.broad.igv.track.Track;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Event signalling need to repaint a track or collection of tracks with existing data.  This event is used when no
 * change in data is indicated, for example a change in track scale or color.
 */
public class RepaintEvent {

    private Collection<Track> tracks;

    public RepaintEvent(Collection<Track> tracks) {
        this.tracks = tracks;
    }

    public RepaintEvent(Track track) {
        this.tracks = Arrays.asList(track);
    }

    public Collection<Track> getTracks() {
        return tracks;
    }
}
