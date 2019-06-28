package org.broad.igv.bedpe;


import htsjdk.samtools.util.Locatable;

import java.awt.*;

public interface BedPE extends Locatable {

    BedPEFeature get();

    double getScore();

    boolean isSameChr();

    void setRow(int row);

    int getRow();

    Color getColor();

    int getThickness();

    void setShape(BedPEShape s);

    BedPEShape getShape();

    String getValueString();

    double getCenterDistance();
}
