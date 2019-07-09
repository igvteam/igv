package org.broad.igv.bedpe;

import java.awt.*;

public class BedPEInterFeature  implements BedPE {


    private final BedPEFeature wrappedFeature;
    private final int order;
    int row;
    BedPEShape shape;

    public BedPEInterFeature(BedPEFeature wrappedFeature, int order) {
        this.wrappedFeature = wrappedFeature;
        this.order = order;
    }

    public BedPEFeature get() {
        return wrappedFeature;
    }

    public String getChr() {
        return this.order == 1 ? wrappedFeature.chr1 : wrappedFeature.chr2;
    }

    public int getStart() {
        return this.order == 1 ? wrappedFeature.start1 : wrappedFeature.start2;
    }

    public int getEnd() {
        return this.order == 1 ? wrappedFeature.end1 : wrappedFeature.end2;
    }

    @Override
    public double getScore() {
        return wrappedFeature.score;
    }

    public boolean isSameChr() {
        return false;
    }

    @Override
    public void setRow(int row) {
        this.row = row;
    }

    @Override
    public int getRow() {
        return row;
    }

    @Override
    public Color getColor() {
        return wrappedFeature.getColor();
    }

    @Override
    public int getThickness() {
        return wrappedFeature.getThickness();
    }

    @Override
    public void setShape(BedPEShape s) {
         this.shape = s;
    }

    @Override
    public BedPEShape getShape() {
        return shape;
    }

    @Override
    public String getValueString() {
        return wrappedFeature.getValueString();
    }

    @Override
    public double getCenterDistance() {
        return 0;
    }

    public String getContig() {
        return getChr();
    }


}
