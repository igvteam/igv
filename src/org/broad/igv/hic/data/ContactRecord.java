package org.broad.igv.hic.data;

/**
 * @author jrobinso
 * @date Aug 3, 2010
 */
public class ContactRecord {

    private int blockNumber;

    /**
     * Bin number in x coordinate
     */
    private int x;

    /**
     * Bin number in y coordinate
     */
    private int y;

    /**
     * Total number of counts
     */
    private int counts;

    public ContactRecord(int block, int x, int y, int counts) {
        this.blockNumber = block;
        this.x = x;
        this.y = y;
        this.counts = counts;
    }

    public void incrementCount() {
        counts++;
    }

    public int getBlockNumber() {
        return blockNumber;
    }

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }

    public int getCounts() {
        return counts;
    }
}
