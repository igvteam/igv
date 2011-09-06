package org.broad.igv.hic;

/**
 * @author jrobinso
 * @date Aug 11, 2010
 */
public class Context {

    private int zoom;
    private int originX = 0;
    private int originY = 0;

    private int visibleWidth;
    private int visibleHeight;
    private int lenX;
    private int lenY;

    public void increment(int dx, int dy) {

       setOrigin(originX + dx, originY + dy);
    }


    public void setOrigin(int x, int y) {
        int maxOriginX = Math.max(0, lenX - visibleWidth);
        int maxOriginY = Math.max(0, lenY - visibleHeight);
        originX = Math.min(maxOriginX, Math.max(0, x));
        originY = Math.min(maxOriginY, Math.max(0, y));
    }

    public int getZoom() {
        return zoom;
    }

    public void setZoom(int zoom) {
        this.zoom = zoom;
    }

    public int getOriginX() {
        return originX;
    }

    public void setOriginX(int originX) {
        this.originX = originX;
    }

    public int getOriginY() {
        return originY;
    }

    public void setOriginY(int originY) {
        this.originY = originY;
    }

    public int getVisibleWidth() {
        return visibleWidth;
    }

    public void setVisibleWidth(int visibleWidth) {
        this.visibleWidth = visibleWidth;
    }

    public int getVisibleHeight() {
        return visibleHeight;
    }

    public void setVisibleHeight(int visibleHeight) {
        this.visibleHeight = visibleHeight;
    }

    public int getLenX() {
        return lenX;
    }

    public void setLenX(int lenX) {
        this.lenX = lenX;
    }

    public int getLenY() {
        return lenY;
    }

    public void setLenY(int lenY) {
        this.lenY = lenY;
    }

}
