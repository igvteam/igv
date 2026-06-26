package org.igv.circview.model;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.UUID;

/**
 * A named set of chords with a shared color and visibility, belonging to a track.
 *
 * <p>Mirrors the chordSet object built by addChords() in circularView.js:
 * {name, trackName, chords, color, trackColor, visible, id}.
 */
public final class ChordSet implements ChordCollection {

    private final String name;
    private final String trackName;
    private final List<Chord> chords;
    private Color color;
    private Color trackColor;
    private boolean visible;
    private final String id;

    public ChordSet(String name, String trackName, List<Chord> chords,
                    Color color, Color trackColor) {
        this.name = name;
        this.trackName = trackName;
        this.chords = new ArrayList<>(chords);
        this.color = color;
        this.trackColor = trackColor;
        this.visible = true;
        this.id = shortId();
    }

    public String getName() {
        return name;
    }

    public String getTrackName() {
        return trackName;
    }

    public List<Chord> getChords() {
        return chords;
    }

    public Color getColor() {
        return color;
    }

    public void setColor(Color color) {
        this.color = color;
    }

    public Color getTrackColor() {
        return trackColor;
    }

    public void setTrackColor(Color trackColor) {
        this.trackColor = trackColor;
    }

    public boolean isVisible() {
        return visible;
    }

    public void setVisible(boolean visible) {
        this.visible = visible;
    }

    public String getId() {
        return id;
    }

    /** Short random id, the Java analogue of guid() in the JS source. */
    static String shortId() {
        return UUID.randomUUID().toString().substring(0, 8);
    }
}
