package org.broad.igv.hic.track;

import org.broad.igv.data.WiggleDataset;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;

import java.awt.*;

/**
 * @author Jim Robinson
 * @date 4/13/12
 */
public class EigenvectorTrack extends AbstractTrack {


    double step;
    double[] data;
    private double dataMax;

    public EigenvectorTrack(String id, String name) {
        super(id, name);
    }

    public void setData(double step, double[] data) {
        this.step = step;
        this.data = data;
        dataMax = 0;
        for (double aData : data) {
            if (Math.abs(aData) > dataMax) dataMax = Math.abs(aData);
        }

    }

    /**
     * Render the track in the supplied rectangle.  It is the responsibility of the track to draw within the
     * bounds of the rectangle.
     *
     * @param context the render context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    public void render(RenderContext context, Rectangle rect) {



        if (data == null) return;

        int h = rect.height / 2;
        Graphics2D g2d = context.getGraphics();
        g2d.setColor(Color.blue.darker());

        int lastXPixel = -1;

        for (int i = 0; i < data.length; i++) {

            int genomicPosition = (int) (step * i);
            int xPixel = context.bpToScreenPixel(genomicPosition);

            if (xPixel > lastXPixel && lastXPixel >= 0) {

                int myh = (int) ((data[i] / dataMax) * h);
                if (data[i] > 0) {
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
}
