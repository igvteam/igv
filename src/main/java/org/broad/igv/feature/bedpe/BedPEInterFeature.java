package org.broad.igv.feature.bedpe;

import java.awt.*;
import java.util.Map;

public class BedPEInterFeature  implements BedPE {


    private final BedPEFeature wrappedFeature;
    private final int order;

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
        return wrappedFeature.isSameChr();
    }

}
