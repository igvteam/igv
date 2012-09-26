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


    double step;
    double[] data;
    private double dataMax;
    private double median;
    HiC hic;
    int currentZoom = -1;

    public EigenvectorTrack(String id, String name, HiC hic) {
          this.hic = hic;
    }

    private void setData(double step, double[] data) {
        this.step = step;
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

    /**
     * Render the track in the supplied rectangle.  It is the responsibility of the track to draw within the
     * bounds of the rectangle.
     *
     * @param g2d     the graphics context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    public void render(Graphics2D g2d, Context context, Rectangle rect) {

        int zoom = hic.zd.getZoom();
        if (zoom != currentZoom) {

            double[] eigen = hic.getEigenvector(0);
            if (eigen == null) return;

            currentZoom = zoom;
            setData(hic.zd.getBinSize(), eigen);
        }

        if (data == null || data.length == 0) return;

        int h = rect.height / 2;
        g2d.setColor(Color.blue.darker());

        int lastXPixel = -1;

        for (int i = 0; i < data.length; i++) {

            if (Double.isNaN(data[i])) continue;

            int genomicPosition = (int) (step * i);

            int xPixel = 0; //context.getScreenPosition (genomicPosition);

            if (xPixel > lastXPixel && lastXPixel >= 0) {

                double x = data[i] - median;
                double max = dataMax - median;

                int myh = (int) ((x / max) * h);
                if (x > 0) {
                    g2d.fillRect(lastXPixel, rect.y + h - myh, (xPixel - lastXPixel), myh);
                } else {
                    g2d.fillRect(lastXPixel, rect.y + h, xPixel - lastXPixel, -myh);
                }
            }
            lastXPixel = xPixel;
        }

    }

    public Renderer getRenderer() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void forceRefresh() {
        currentZoom = -1;
    }
}
