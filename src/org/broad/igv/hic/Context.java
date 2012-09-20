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
    private int genomicOrigin = 0;
    private double scale;

    private int binOrigin = 0;

    public Context(Chromosome chromosome) {
        this.chromosome = chromosome;
    }

    public int getBinCount(int binSize) {
        return (getChrLength() - getGenomicOrigin()) / binSize + 1;
    }

    public int getBinOrigin() {
        return binOrigin;
    }

    public void setBinOrigin(int binOrigin) {
        this.binOrigin = binOrigin;
    }

    public void setGenomicOrigin(int x) {
        genomicOrigin = Math.max(0, x);

    }

    public int getZoom() {
        return zoom;
    }

    public void setZoom(int zoom) {
        this.zoom = zoom;
    }

    public void setZoom(int zoom, double scale) {
        this.scale = scale;
        this.zoom = zoom;
    }

    public int getBinNumber() {
        return (int) (genomicOrigin / scale);
    }

    public int getGenomicOrigin() {
        return genomicOrigin;
    }

    public int getChrLength() {
        return chromosome.getLength();
    }

    public double getScale() {
        return scale;
    }

    /**
     * Return the screen position corresponding to the chromosomal position.
     *
     * @param chromosomePosition
     * @return
     */
    public int getScreenPosition(double chromosomePosition) {
        return (int) ((chromosomePosition - genomicOrigin) / scale);
    }

    public double getChromosomePosition(int screenPosition) {
        return genomicOrigin + screenPosition * scale;
    }

    public Chromosome getChromosome() {
        return chromosome;
    }

}
