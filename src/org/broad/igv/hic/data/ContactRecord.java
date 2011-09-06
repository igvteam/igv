package org.broad.igv.hic.data;

/**
 * @author jrobinso
 * @date Aug 3, 2010
 */
public class ContactRecord {

    private int blockNumber;
    private int x;
    private int y;
    private short counts;

    public ContactRecord(int block, int x, int bin2, short counts) {
        this.blockNumber = block;
        this.x = x;
        this.y = bin2;
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

    public short getCounts() {
        return counts;
    }
}
