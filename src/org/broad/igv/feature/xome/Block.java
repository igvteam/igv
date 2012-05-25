package org.broad.igv.feature.xome;


import org.broad.tribble.Feature;

/**
 * @author Jim Robinson
 * @date 5/24/12
 */
public class Block implements Feature {
    int idx;
    private int genomeStart;
    private int genomeEnd;
    private int exomeStart;

    public Block(int idx, int start, int end, int xomeStart) {
        this.idx = idx;
        this.genomeStart = start;
        this.genomeEnd = end;
        this.exomeStart = xomeStart;
    }

    public int getIdx() {
        return idx;
    }

    public int getGenomeStart() {
        return genomeStart;
    }

    public int getGenomeEnd() {
        return genomeEnd;
    }

    public int getExomeStart() {
        return exomeStart;
    }

    public int getExomeEnd() {
        return exomeStart + (genomeEnd - genomeStart);
    }

    public void extend(int x) {
        if (x > getGenomeEnd()) genomeEnd = x;
    }

    public int getLength() {
        return getGenomeEnd() - getGenomeStart();
    }


    public String toString() {
        return "Block " + idx + " [" + genomeStart + ", " + genomeEnd + ", " + exomeStart + "]";
    }

    @Override
    public String getChr() {

        return null;
    }

    @Override
    public int getStart() {
        return exomeStart;
    }

    @Override
    public int getEnd() {
        return exomeStart + (genomeEnd - genomeStart);
    }
}
