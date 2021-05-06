package org.broad.igv.sam;

public class BaseModifications {

    public char base;
    public char strand;
    public String modification;
    public int[] positions;

    public BaseModifications(char base, char strand, String modification, int[] positions) {
        this.base = base;
        this.strand = strand;
        this.modification = modification;
        this.positions = positions;
    }

    public static class Mod {
        char base;
        char strand;
        byte likelihood;
        public Mod(char base, char strand, byte likelihood) {
            this.base = base;
            this.strand = strand;
            this.likelihood = likelihood;
        }
    }
}
