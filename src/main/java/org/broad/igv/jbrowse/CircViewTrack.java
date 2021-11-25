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

    public String toJson() {
        StringBuffer buf = new StringBuffer();
        buf.append("{");
        buf.append(JsonUtils.toJson("name", name));
        buf.append(",");
        buf.append(JsonUtils.toJson("color", color));
        buf.append(",\"chords\":");
        buf.append("[");
        boolean first = true;
        for(Chord c : chords) {
            if(!first) {
                buf.append(",");
            }
            buf.append(c.toJson());
            first = false;
        }
        buf.append("]");
        buf.append("}");
        return buf.toString();
    }
}
