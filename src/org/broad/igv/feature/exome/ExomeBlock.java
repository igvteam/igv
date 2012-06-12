package org.broad.igv.feature.exome;


import org.broad.tribble.Feature;

/**
 * @author Jim Robinson
 * @date 5/24/12
 */
public class ExomeBlock implements Feature {
    int idx;
    private int genomeStart;
    private int exomeStart;
    private int length;
    private int leftPixel;
    private int rightPixel;

    public ExomeBlock(int idx, int genomeStart, int xomeStart, int length) {
        this.idx = idx;
        this.genomeStart = genomeStart;
        this.exomeStart = xomeStart;
        this.length = length;
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
        return genomeStart + length;
    }

    public int getExomeStart() {
        return exomeStart;
    }

    public int getExomeEnd() {
        return exomeStart + length;
    }

    public void extend(int x) {
        if (x > genomeStart + length) {
            length = x - genomeStart;
        };
    }

    public int getLength() {
        return getGenomeEnd() - getGenomeStart();
    }


    public String toString() {
        return "Block " + idx + " [" + genomeStart + ", " + getGenomeEnd() + ", " + exomeStart + ", " + length +"]";
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
        return genomeStart + length;
    }

    public int compareGenomePosition(double genomicPosition) {
        if(genomicPosition < genomeStart) {
            return -1;
        }
        else if(genomicPosition >= genomeStart + length) {
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
        else if(exomePosition >= (exomeStart + length)) {
            return 1;
        }
        else {
            return 0;
        }
    }
}
