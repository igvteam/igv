package org.broad.igv.ucsc;

import java.util.List;

public class TrackConfigGroup {

    public String label;
    public List<TrackConfig> tracks;

    public TrackConfigGroup(String label, List<TrackConfig> tracks) {
        this.label = label;
        this.tracks = tracks;
    }
}
