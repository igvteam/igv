package org.broad.igv.jbrowse;

class Mate {
    String refName;
    int start;
    int end;

    public Mate(String refName, int start, int end) {
        this.refName = refName;
        this.start = start;
        this.end = end;
    }
}
