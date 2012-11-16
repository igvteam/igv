package org.broad.igv.hic.data;

/**
 * @author jrobinso
 * @date Aug 3, 2010
 */
public class ContactRecord implements Comparable<ContactRecord> {

    /**
     * Bin number in x coordinate
     */
    private int binX;

    /**
     * Bin number in y coordinate
     */
    private int binY;

    /**
     * Total number of counts, or cumulative score
     */
    private float counts;

    public ContactRecord(int binX, int binY, float counts) {
        this.binX = binX;
        this.binY = binY;
        this.counts = counts;
    }

    public void incrementCount(float score) {
        counts += score;
    }


    public int getBinX() {
        return binX;
    }

    public int getBinY() {
        return binY;
    }

    public float getCounts() {
        return counts;
    }

    @Override
    public int compareTo(ContactRecord contactRecord) {
        if(this.binX != contactRecord.binX) {
            return binX - contactRecord.binX;
        }
        else if(this.binY != contactRecord.binY) {
            return binY - contactRecord.binY;
        }
        else return 0;
    }
}
