package org.igv.feature.basepair;

import htsjdk.tribble.Feature;

import java.awt.*;

/**
 * Created by jrobinson on 3/1/16.
 */
public class BasePairFeature implements Feature{

    String chr;
    int startLeft;
    int startRight;
    int endLeft;
    int endRight;
    int colorIndex;

    public BasePairFeature(String chr,  int startLeft, int startRight, int endLeft, int endRight,  int colorIndex) {
        this.chr = chr;
        this.colorIndex = colorIndex;
        this.endLeft = endLeft;
        this.endRight = endRight;
        this.startLeft = startLeft;
        this.startRight = startRight;
    }

    @Override
    public String getChr() {
        return chr;
    }

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public int getStart() {
        return startLeft;
    }

    @Override
    public int getEnd() {
        return endRight;
    }

    public int getStartLeft() { return startLeft; }

    public int getStartRight() { return startRight; }

    public int getEndLeft() { return endLeft; }

    public int getEndRight() { return endRight; }

    public int getColorIndex() { return colorIndex; }

    public String toString() {
        return getChr() + "\t" + startLeft + "\t" + startRight + "\t" + endLeft + "\t" + endRight;
    }
}
