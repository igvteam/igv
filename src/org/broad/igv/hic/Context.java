package org.broad.igv.hic;

//import org.broad.igv.hic.data.Chromosome;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.track.HiCFixedGridAxis;
import org.broad.igv.hic.track.HiCGridAxis;

/**
 * @author jrobinso
 * @date Aug 11, 2010
 */
public class Context {

    private Chromosome chromosome;
    private int zoom = 4;

    private double scaleFactor = 1;

    private int binOrigin = 0;

    public Context(Chromosome chromosome) {
        this.chromosome = chromosome;
    }

    public int getBinOrigin() {
        return binOrigin;
    }

    public void setBinOrigin(int binOrigin) {
        this.binOrigin = binOrigin;
    }

    public int getZoom() {
        return zoom;
    }

    public void setZoom(int zoom) {
        this.zoom = zoom;
    }

    public void setZoom(int zoom, double scale) {
        this.scaleFactor = scale;
        this.zoom = zoom;
    }

    public int getChrLength() {
        return chromosome.getLength();
    }


    public Chromosome getChromosome() {
        return chromosome;
    }

    public void setScaleFactor(double scale) {
        this.scaleFactor = scale;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }
}
