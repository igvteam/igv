package org.broad.igv.hic.track;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.data.WiggleDataset;
import org.broad.igv.hic.Context;
import org.broad.igv.hic.HiC;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.util.collections.DoubleArrayList;

import java.awt.*;

/**
 * @author Jim Robinson
 * @date 4/13/12
 */
public class EigenvectorTrack extends HiCTrack {


    public static final Color COLOR = Color.blue.darker();
    double scale;
    double[] data;
    private double dataMax;
    private double median;
    HiC hic;
    int currentZoom = -1;

    public EigenvectorTrack(String id, String name, HiC hic) {
        this.hic = hic;
    }

    private void setData(double scale, double[] data) {
        this.scale = scale;
        this.data = data;

        DoubleArrayList tmp = new DoubleArrayList(data.length);
        for (int i = 0; i < data.length; i++) {
            if (!Double.isNaN(data[i])) {
                tmp.add(data[i]);
            }
        }
        double[] tmpArray = tmp.toArray();
        this.median = StatUtils.percentile(tmpArray, 50);
        dataMax = 0;
        for (double aData : tmpArray) {
            if (Math.abs(aData) > dataMax) dataMax = Math.abs(aData);
        }
    }

    public Color getColor() {
        return COLOR;
    }

    /**
     * Render the track in the supplied rectangle.  It is the responsibility of the track to draw within the
     * bounds of the rectangle.
     *
     * @param g2d  the graphics context
     * @param rect the track bounds, relative to the enclosing DataPanel bounds.
     */
    public void render(Graphics2D g2d, Context context, Rectangle rect) {

        int zoom = hic.zd.getZoom();
        if (zoom != currentZoom) {

            double[] eigen = hic.getEigenvector(0);
            if (eigen == null) return;

            currentZoom = zoom;
            setData(context.getScaleFactor(), eigen);
        }

        if (data == null || data.length == 0) return;

        int h = rect.height / 2;
        g2d.setColor(COLOR);

        for (int bin = context.getBinOrigin(); bin < data.length; bin++) {

            if (Double.isNaN(data[bin])) continue;

            int xPixelLeft = rect.x + (int) ((bin - context.getBinOrigin()) * scale); //context.getScreenPosition (genomicPosition);
            int xPixelRight = rect.x + (int) ((bin + 1 - context.getBinOrigin()) * scale);


            if (xPixelRight < rect.x) {
                continue;
            } else if (xPixelLeft > rect.x + rect.width) {
                break;
            }

            double x = data[bin] - median;
            double max = dataMax - median;

            int myh = (int) ((x / max) * h);
            if (x > 0) {
                g2d.fillRect(xPixelLeft, rect.y + h - myh, (xPixelRight - xPixelLeft), myh);
            } else {
                g2d.fillRect(xPixelLeft, rect.y + h, (xPixelRight - xPixelLeft), -myh);
            }

        }

    }

    public String getName() {
        return "eigenvector";
    }

    public Renderer getRenderer() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void forceRefresh() {
        currentZoom = -1;
    }
}
