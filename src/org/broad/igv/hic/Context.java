package org.broad.igv.hic;

/**
 * @author jrobinso
 * @date Aug 11, 2010
 */
public class Context {

    private int zoom = 4;
    private int origin = 0;

    private int visibleWidth;

    /**
     * Total length of the chromosome
     */
    private int chrLength;
    private double scale;


    public void increment(int origin) {
       setOrigin(origin);
    }


    public void setOrigin(int x) {
        int maxOrigin = Math.max(0, chrLength - visibleWidth);
        origin = Math.min(maxOrigin, Math.max(0, x));
    }

    public int getZoom() {
        return zoom;
    }

    public void setZoom(int zoom, double v) {
        this.scale = v;
        this.zoom = zoom;
    }

    public int getOrigin() {
        return origin;
    }

    public int getVisibleWidth() {
        return visibleWidth;
    }

    public void setVisibleWidth(int visibleWidth) {
        this.visibleWidth = visibleWidth;
    }

    public int getChrLength() {
        return chrLength;
    }

    public void setChrLength(int lenX) {
        this.chrLength = lenX;
    }

    public double getScale() {
        return scale;  //To change body of created methods use File | Settings | File Templates.
    }

    /**
     * Return the screen position corresponding to the chromosomal position.
     *
     * @param chromosomePosition
     * @return
     */
    public int getScreenPosition(double chromosomePosition) {
        return (int) ((chromosomePosition - origin) / getScale());
    }
}
