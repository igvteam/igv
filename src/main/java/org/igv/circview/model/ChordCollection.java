package org.igv.circview.model;

import java.awt.Color;
import java.util.List;

/**
 * A named, colorable, hideable collection of chords. Implemented by both
 * {@link ChordSet} (the flat view) and {@link Track} (the grouped view) so the
 * controls can drive either without caring which is active.
 */
public interface ChordCollection {

    String getName();

    Color getColor();

    void setColor(Color color);

    boolean isVisible();

    void setVisible(boolean visible);

    /** All chords in this collection. */
    List<Chord> getChords();
}
