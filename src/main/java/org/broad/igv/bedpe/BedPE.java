package org.broad.igv.bedpe;

import org.broad.igv.feature.Locatable;

public interface BedPE extends Locatable {

    BedPEFeature get();

    double getScore();

    boolean isSameChr();

    void setRow(int row);

    int getRow();
}
