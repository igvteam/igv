package org.igv.circview.model;

/**
 * The far end of a chord: a genomic region {refName, start, end}.
 * Mirrors the "mate" object on a chord feature in circularView.js.
 */
public final class Mate {

    private final String refName;
    private final long start;
    private final long end;

    public Mate(String refName, long start, long end) {
        this.refName = refName;
        this.start = start;
        this.end = end;
    }

    public String getRefName() {
        return refName;
    }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return end;
    }
}
