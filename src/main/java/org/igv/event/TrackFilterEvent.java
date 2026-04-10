package org.igv.event;

import org.igv.sample.SampleFilter;

public final class TrackFilterEvent implements IGVEvent {

    private SampleFilter sampleFilter;

    public TrackFilterEvent(SampleFilter sampleFilter) {
        this.sampleFilter = sampleFilter;
    }

    public SampleFilter getFilter() {
        return sampleFilter;
    }
}
