package org.broad.igv.hic;

//import org.broad.igv.hic.data.Chromosome;
import org.broad.igv.feature.Chromosome;

/**
 * @author jrobinso
 * @date Aug 11, 2010
 */
public class Context {

    private Chromosome chromosome;
    private int zoom = 4;
    private int origin = 0;
    private double scale;

    public Context(Chromosome chromosome) {
        this.chromosome = chromosome;
    }


    public void setOrigin(int x) {
        origin = Math.max(0, x);

    }

    public int getZoom() {
        return zoom;
    }

    public void setZoom(int zoom, double scale) {
        this.scale = scale;
        this.zoom = zoom;
    }

    public int getOrigin() {
        return origin;
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
        return (int) ((chromosomePosition - origin) / scale);
    }

    public double getChromosomePosition(int screenPosition) {
        return origin + screenPosition * scale;
    }

    public Chromosome getChromosome() {
        return chromosome;
    }

}
