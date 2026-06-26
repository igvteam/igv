package org.igv.circview.ui;

import org.igv.circview.model.Chord;

/**
 * Callback invoked when a chord is clicked. The Java analogue of the
 * onChordClick config option in circularView.js.
 */
@FunctionalInterface
public interface ChordClickListener {
    void onChordClick(Chord feature);
}
