package org.igv.feature.tribble;

import org.igv.track.TrackProperties;
import org.igv.track.TrackType;

import java.util.*;

//  TODO -- why do we need this class and TrackProperties?

public class FeatureFileHeader {

    /* An object to collection track properties, if specified in the feature file. */
    private TrackProperties trackProperties;

    private TrackType trackType;

    public FeatureFileHeader() {
    }

    public FeatureFileHeader(TrackProperties trackProperties) {
        this.trackProperties = trackProperties;
    }

    public void setTrackProperties(TrackProperties trackProperties) {
        this.trackProperties = trackProperties;
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public void setTrackType(TrackType trackType) {
        this.trackType = trackType;
    }

    public TrackType getTrackType() {
        return trackType;
    }

}
