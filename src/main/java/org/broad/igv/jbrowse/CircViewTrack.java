package org.broad.igv.jbrowse;

class CircViewTrack {
    String name;
    String color;
    Chord[] chords;

    public CircViewTrack(Chord[] chords, String name, String color) {
        this.name = name;
        this.color = color;
        this.chords = chords;
    }
}
