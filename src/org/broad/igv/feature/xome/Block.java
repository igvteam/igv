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
    private int leftPixel;
    private int rightPixel;

    public Block(int idx, int start, int end, int xomeStart) {
        this.idx = idx;
        this.genomeStart = start;
        this.genomeEnd = end;
        this.exomeStart = xomeStart;
    }

    public void setScreenBounds(int leftPixel, int rightPixel) {
        this.leftPixel = leftPixel;
        this.rightPixel = rightPixel;
    }

    public int getLeftPixel() {
        return leftPixel;
    }

    public int getRightPixel() {
        return rightPixel;
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
        if (x > genomeEnd) genomeEnd = x;
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
        return genomeStart;
    }

    @Override
    public int getEnd() {
        return genomeEnd;
    }

    public int compareGenomePosition(double genomicPosition) {
        if(genomicPosition < genomeStart) {
            return -1;
        }
        else if(genomicPosition >= genomeEnd) {
            return 1;
        }
        else {
            return 0;
        }

    }

    public int compareExomePosition(int exomePosition) {
        if(exomePosition < exomeStart) {
            return -1;
        }
        else if(exomePosition >= (exomeStart + (genomeEnd - genomeStart))) {
            return 1;
        }
        else {
            return 0;
        }
    }
}
