package org.broad.igv.ui.panel;

import org.broad.igv.track.Track;

import java.util.ArrayList;
import java.util.List;

public class TrackWrapper<T extends Track> {
    private T track;

    public TrackWrapper(T track) {
        this.track = track;
    }

    public String toString() {
        return track.getName();
    }

    public T getTrack() {
        return this.track;
    }

    /**
     *
     * @param tracks
     * @return
     */
    public static List<TrackWrapper> wrapTracks(Iterable<? extends Track> tracks) {
        ArrayList<TrackWrapper> wrappers = new ArrayList<TrackWrapper>();
        for (Track t : tracks) {
            TrackWrapper trackWrapper = new TrackWrapper(t);
            wrappers.add(trackWrapper);
        }
        return wrappers;
    }

}