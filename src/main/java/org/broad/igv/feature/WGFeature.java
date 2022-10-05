package org.broad.igv.feature;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.WindowFunction;

import java.awt.*;

public class WGFeature implements IGVFeature {

    IGVFeature wrappedFeature;
    int start;
    int end;

    public WGFeature(IGVFeature wrappedFeature, int start, int end) {
        this.wrappedFeature = wrappedFeature;
        this.start = start;
        this.end = end;
    }


    @Override
    public Color getColor() {
        return wrappedFeature.getColor();
    }

    @Override
    public float getScore() {
        return wrappedFeature.getScore();
    }

    @Override
    public String getName() {
        return wrappedFeature.getName();
    }

    @Override
    public String getContig() {
        return Globals.CHR_ALL;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return end;
    }

    @Override
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        // TODDO -- translate position to chr position
        return wrappedFeature.getValueString(0, mouseX, windowFunction);
    }
}
