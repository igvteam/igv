package org.broad.igv.sam;

/**
 * Created by jrobinso on 3/31/16.
 */
public class SpliceGap extends Gap {

    int flankingLeft;
    int flankingRight;

    public SpliceGap(int start, int nBases, char type, int flankingLeft, int flankingRight) {
        super(start, nBases, type);
        this.flankingLeft = flankingLeft;
        this.flankingRight = flankingRight;
    }

    public int getFlankingLeft() {
        return flankingLeft;
    }

    public int getFlankingRight() {
        return flankingRight;
    }
}
