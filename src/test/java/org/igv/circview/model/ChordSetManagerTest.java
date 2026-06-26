package org.igv.circview.model;

import org.junit.Test;

import java.awt.Color;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

public class ChordSetManagerTest {

    private static Chord chord(String ref) {
        return new Chord(ref + "-id", ref, 1, 2, new Mate(ref, 3, 4), null);
    }

    private static ChordSet set(String name, String trackName, Chord... chords) {
        return new ChordSet(name, trackName, List.of(chords), Color.BLACK, Color.BLACK);
    }

    @Test
    public void replacesChordSetWithSameName() {
        ChordSetManager m = new ChordSetManager();
        m.addChordSet(set("A", "trackA", chord("1")));
        m.addChordSet(set("A", "trackA", chord("2"), chord("3")));

        assertEquals(1, m.getChordSets().size());
        assertEquals(2, m.getChordSet("A").getChords().size());
    }

    @Test
    public void groupsChordSetsByTrack() {
        ChordSetManager m = new ChordSetManager();
        // Track name is derived externally; here both sets share "trackA".
        m.addChordSet(set("region1", "trackA", chord("1")));
        m.addChordSet(set("region2", "trackA", chord("2")));

        assertEquals(2, m.getChordSets().size());
        assertEquals(1, m.getTracks().size());

        Track track = m.getTrack("trackA");
        assertNotNull(track);
        assertEquals(2, track.getChordSets().size());
        // Track.chords() concatenates across its chord sets.
        assertEquals(2, track.getChords().size());
    }

    @Test
    public void separateTracksStaySeparate() {
        ChordSetManager m = new ChordSetManager();
        m.addChordSet(set("region1", "trackA", chord("1")));
        m.addChordSet(set("region2", "trackB", chord("2")));

        assertEquals(2, m.getTracks().size());
        assertEquals(1, m.getTrack("trackA").getChords().size());
        assertEquals(1, m.getTrack("trackB").getChords().size());
    }

    @Test
    public void clearChordsEmptiesEverything() {
        ChordSetManager m = new ChordSetManager();
        m.addChordSet(set("A", "trackA", chord("1")));
        m.clearChords();
        assertEquals(0, m.getChordSets().size());
        assertEquals(0, m.getTracks().size());
        assertNull(m.getChordSet("A"));
    }
}
