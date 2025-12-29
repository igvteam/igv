package org.broad.igv.track;

import org.broad.igv.feature.Locus;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.List;

/**
 * @author jrobinso
 * @date Sep 28, 2010
 */
public class LoadedDataInterval<T> {


    public Locus range;
    private T scores;
    int zoom = -1;

    public LoadedDataInterval(String chr, int start, int end, T scores) {

        range = new Locus(chr, start, end);
        this.scores = scores;
    }

    public LoadedDataInterval(String chr, int start, int end, int zoom, T scores) {

        range = new Locus(chr, start, end);
        this.zoom = zoom;
        this.scores = scores;
    }
    public boolean contains(String chr, int start, int end) {
        return range.contains(chr, start, end);
    }

    public boolean contains(String chr, int start, int end, int zoom) {
        return (this.zoom == -1 || this.zoom == zoom) && range.contains(chr, start, end);
    }

    public boolean contains(ReferenceFrame frame) {
        String chr = frame.getChrName();
        int start = (int) frame.getOrigin();
        int end = (int) frame.getEnd();
        int zoom = frame.getZoom();
        return this.contains(chr, start, end, zoom);
    }

    public T getFeatures() {
        return scores;
    }
}
