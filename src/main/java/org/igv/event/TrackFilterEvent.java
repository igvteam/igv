package org.igv.event;

import org.igv.util.TrackFilter;

public final class TrackFilterEvent implements IGVEvent {

    private TrackFilter trackFilter;

    public TrackFilterEvent(TrackFilter trackFilter) {
        this.trackFilter = trackFilter;
    }

    public TrackFilter getFilter() {
        return trackFilter;
    }
}
