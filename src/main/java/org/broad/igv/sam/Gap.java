package org.broad.igv.sam;

/**
 * Created by jrobinso on 3/22/16.
 */
public class Gap {
    int start;
    int nBases;
    char type;

    public Gap(int start, int nBases, char type) {
        this.start = start;
        this.nBases = nBases;
        this.type = type;
    }

    public int getStart() {
        return start;
    }

    public int getnBases() {
        return nBases;
    }

    public char getType() {
        return type;
    }
}
