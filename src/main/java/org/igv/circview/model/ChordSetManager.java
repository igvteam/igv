package org.igv.circview.model;

import java.util.ArrayList;
import java.util.List;

/**
 * Maintains the set of chords as both a flat list of {@link ChordSet}s and a
 * grouped list of {@link Track}s. The view chooses which list to render based on
 * its "group by track" flag.
 *
 * <p>Direct port of chordSetManager.js.
 */
public final class ChordSetManager {

    private final List<Track> tracks = new ArrayList<>();
    private final List<ChordSet> chordSets = new ArrayList<>();

    public void addChordSet(ChordSet chordSet) {
        // If a chord set with this name exists, replace it (same track, same region).
        chordSets.removeIf(g -> g.getName().equals(chordSet.getName()));
        chordSets.add(chordSet);

        Track track = null;
        for (Track t : tracks) {
            if (chordSet.getTrackName().equals(t.getName())) {
                track = t;
                break;
            }
        }
        if (track != null) {
            track.getChordSets().removeIf(cs -> cs.getName().equals(chordSet.getName()));
            track.getChordSets().add(chordSet);
        } else {
            tracks.add(new Track(chordSet));
        }
    }

    public void clearChords() {
        tracks.clear();
        chordSets.clear();
    }

    public List<Track> getTracks() {
        return tracks;
    }

    public List<ChordSet> getChordSets() {
        return chordSets;
    }

    public Track getTrack(String name) {
        for (Track t : tracks) {
            if (t.getName().equals(name)) {
                return t;
            }
        }
        return null;
    }

    public ChordSet getChordSet(String name) {
        for (ChordSet cs : chordSets) {
            if (cs.getName().equals(name)) {
                return cs;
            }
        }
        return null;
    }
}
