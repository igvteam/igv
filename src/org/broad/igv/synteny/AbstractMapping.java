package org.broad.igv.synteny;

public abstract class AbstractMapping implements Mapping {
    private String name;
    protected String fromChr;
    protected int fromStart;
    protected int fromEnd;
    protected String toChr;
    protected int toStart;
    protected int toEnd;
    protected boolean direction;
    double scaleFactor;

    public void setParameters(String name, String fromChr, int fromStart, int fromEnd, String fromDir,
                              String toChr, int toStart, int toEnd, String toDir) {
        this.name = name;
        this.fromChr = fromChr;
        this.fromStart = fromStart;
        this.fromEnd = fromEnd;
        this.toChr = toChr;
        this.toStart = toStart;
        this.toEnd = toEnd;
        this.direction = fromDir.equals(toDir);
        this.scaleFactor = ((double) (toEnd - toStart) / (fromEnd - fromStart));
    }

    public String toString() {
        return this.name + " " + this.fromChr + ":" + this.fromStart + "-" + this.fromEnd + " -> " + this.toChr + ":" + this.toStart + "-" + this.toEnd;
    }

    public boolean containsFromPosition(int fromPosition) {
        return (fromPosition >= this.fromStart) && (fromPosition <= this.fromEnd);
    }

    public String getName() {
        return this.name;
    }

    public String getFromChr() {
        return this.fromChr;
    }

    public int getFromStart() {
        return this.fromStart;
    }

    public int getFromEnd() {
        return this.fromEnd;
    }

    public String getToChr() {
        return this.toChr;
    }

    public int getToStart() {
        return this.toStart;
    }

    public int getToEnd() {
        return this.toEnd;
    }

    public boolean getDirection() {
        return this.direction;
    }
}