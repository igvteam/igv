package org.igv.feature.tribble;

import org.igv.track.TrackProperties;
import org.igv.track.DataType;

//  TODO -- why do we need this class and TrackProperties?

public class FeatureFileHeader {

    /* An object to collection track properties, if specified in the feature file. */
    private TrackProperties trackProperties;

    private DataType dataType;

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

    public void setTrackType(DataType dataType) {
        this.dataType = dataType;
    }

    public DataType getTrackType() {
        return dataType;
    }

}
