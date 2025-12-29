/*
 * TrackSet.java
 *   A TrackSet is a collection of tracks that share a common location axis.
 */
package org.igv.track;

import java.util.List;

/**
 * @author jrobinso
 */
public class TrackSet {

    private String name;
    List<Track> dataTracks;

    public TrackSet(List<Track> dataTracks) {
        this.dataTracks = dataTracks;
    }

    public List<Track> getTracks() {
        return dataTracks;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public boolean isEmpty() {
        return dataTracks == null || dataTracks.isEmpty();
    }
}


