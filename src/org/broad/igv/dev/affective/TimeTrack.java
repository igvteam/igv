package org.broad.igv.dev.affective;

import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;

import java.awt.*;

/**
 * @author Jim Robinson
 * @date 1/11/12
 */
public class TimeTrack extends AbstractTrack {

    double timeStep = 1.0 / 8;

    public TimeTrack(String id) {
        super(id);
    }


    @Override
    public int getHeight() {
        return 20;
    }

    /**
     * Render the track in the supplied rectangle.  It is the responsibility of the track to draw within the
     * bounds of the rectangle.
     *
     * @param context the render context
     * @param rect    the track bounds, relative to the enclosing DataPanel bounds.
     */
    public void render(RenderContext context, Rectangle rect) {

        double start = context.getReferenceFrame().getOrigin();
        double end = context.getReferenceFrame().getEnd();
        double seconds = (start - end) * timeStep;

        // Determine step sizes
        double secsPerPixel = context.getScale() * timeStep;
        double minsPerPixel = secsPerPixel / 60;
        double hoursPerPixel = minsPerPixel / 60;

        // Mark ticks
        Graphics2D g = context.getGraphic2DForColor(Color.black);

        // start and end time in hours, rounded
        int startHour = (int) ((start * timeStep) / 3600);
        int endHour = (int) ((end * timeStep) / 3600) + 1;

        double originHour = (start * timeStep) / 3600;
        System.out.println(startHour + " - " + endHour);
        for (double h = startHour; h < endHour; h++) {
            double pixel = (int) ((h - originHour) / hoursPerPixel);
            //if(pixel < 0) {
            //    System.out.println(h + " " + start + " " + (h - start ) + " " + pixel);
            //}
            if (pixel > rect.x + rect.width) {
                break;
            }
            if (pixel >= 0) {
                int base = rect.y + rect.height - 15;
                //GraphicUtils.drawCenteredChar(String.valueOf((int) h),
                g.drawLine((int) pixel, rect.y, (int) pixel, rect.y + rect.height);

                // If room for 1/4 hours
                if (15 / minsPerPixel > 10) {
                      for(int mm=0; mm < 60; mm+= 15) {
                          double dx = mm / minsPerPixel;
                          double mPixel = pixel + dx;
                          g.drawLine((int) mPixel, rect.y + rect.height/2, (int) mPixel, rect.y + rect.height);

                      }
                }

                // If room for minutes do minutes
                if (minsPerPixel < 0.5) {
                      for(int m=1; m<60; m++) {
                          double dx = m / minsPerPixel;
                          double mPixel = pixel + dx;
                          g.drawLine((int) mPixel, rect.y + 3*rect.height/4, (int) mPixel, rect.y + rect.height);

                      }
                }

            }
        }


    }

    public Renderer getRenderer() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
