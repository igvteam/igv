package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.TrackConfig;

import java.util.ArrayList;
import java.util.List;

public class TrackConfigGroup {

    public int priority;
    public String name;
    public String label;
    public boolean defaultOpen;
    public List<TrackConfig> tracks;

    public TrackConfigGroup(String name, String label, int priority, boolean defaultOpen) {
        this.name = name;
        this.priority = priority;
        this.label = label;
        this.defaultOpen = defaultOpen;
        this.tracks = new ArrayList<>();
    }

    public boolean isEmpty() {
        return tracks.isEmpty();
    }
}
