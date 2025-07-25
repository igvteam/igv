package org.broad.igv.event;

import org.broad.igv.util.Filter;

public final class TrackFilterEvent implements IGVEvent {

    private Filter filter;

    public TrackFilterEvent(Filter filter) {
        this.filter = filter;
    }

    public Filter getFilter() {
        return filter;
    }
}
