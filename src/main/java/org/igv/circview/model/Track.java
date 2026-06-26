package org.igv.circview.model;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

/**
 * A track groups one or more {@link ChordSet}s under a common name and color.
 * Port of IGVTrack in chordSetManager.js.
 */
public final class Track implements ChordCollection {

    private final String name;
    private Color color;
    private boolean visible;
    private final List<ChordSet> chordSets;
    private final String id;

    public Track(ChordSet chordSet) {
        this.name = chordSet.getTrackName();
        this.color = chordSet.getTrackColor();
        this.visible = true;
        this.chordSets = new ArrayList<>();
        this.chordSets.add(chordSet);
        this.id = ChordSet.shortId();
    }

    public String getName() {
        return name;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
    }

    public boolean isVisible() {
        return visible;
    }

    public void setVisible(boolean visible) {
        this.visible = visible;
    }

    public List<ChordSet> getChordSets() {
        return chordSets;
    }

    public String getId() {
        return id;
    }

    /**
     * All chords across this track's chord sets.
     * Mirrors the {@code get chords()} accessor on IGVTrack.
     */
    @Override
    public List<Chord> getChords() {
        if (chordSets.size() == 1) {
            return chordSets.get(0).getChords();
        }
        List<Chord> all = new ArrayList<>();
        for (ChordSet cs : chordSets) {
            all.addAll(cs.getChords());
        }
        return all;
    }
}
