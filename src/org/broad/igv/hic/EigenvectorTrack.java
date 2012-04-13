package org.broad.igv.hic;

import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;

import java.awt.*;

/**
 * @author Jim Robinson
 * @date 4/13/12
 */
public class EigenvectorTrack extends AbstractTrack {

    double [] data;

    public EigenvectorTrack(String id, String name) {
        super(id, name);
    }

    public void setData(double[] data) {
        this.data = data;
    }

    /**
     * Render the track in the supplied rectangle.  It is the responsibility of the track to draw within the
     * bounds of the rectangle.
     *
     * @param context the render context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    public void render(RenderContext context, Rectangle rect) {
        int h = rect.height / 2;
        Graphics2D g2d = context.getGraphic2DForColor(Color.blue.darker());
        g2d.setColor(Color.blue.darker());
        double data_max = 0;
        for (double aData : data) {
            if (Math.abs(aData) > data_max) data_max = Math.abs(aData);
        }
        for (int i = 0; i < data.length; i++) {
            int myh = (int) ((data[i] / data_max) * h);
            if (data[i] > 0) {
                g2d.fillRect(i, h - myh, 1, myh);
            } else {
                System.out.println(h + " " + myh);
                g2d.fillRect(i, h, 1, -myh);
            }
        }

    }

    public Renderer getRenderer() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
