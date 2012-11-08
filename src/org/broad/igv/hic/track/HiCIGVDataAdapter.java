package org.broad.igv.hic.track;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.hic.HiC;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.DataTrack;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 11/8/12
 *         Time: 10:16 AM
 */
public class HiCIGVDataAdapter extends HiCDataAdapter {

    DataTrack igvTrack;


    public HiCIGVDataAdapter(HiC hic, DataTrack igvTrack) {
        super(hic);
        this.igvTrack = igvTrack;
    }



    public double getMax() {
        return igvTrack.getDataRange().getMaximum();
    }

    public String getName() {
        return igvTrack.getName();
    }

    public Color getColor() {
        return igvTrack.getColor();
    }

    public boolean isLogScale() {
        return igvTrack.getDataRange().isLog();
    }

    public Color getAltColor() {
        return igvTrack.getAltColor();
    }

    public DataRange getDataRange() {
        return igvTrack.getDataRange();
    }

    protected List<LocusScore> getLocusScores(String chr, int zoom, int gStart, int gEnd) {
        return igvTrack.getSummaryScores("chr" + chr, gStart, gEnd, zoom);
    }
}
