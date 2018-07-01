package org.broad.igv.feature.bedpe;

import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.util.List;


/**
 * Created by jrobinso on 6/29/18.
 */
public class PEArcRenderer {


    public void render(List<BedPEFeature> featureList, RenderContext context, Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        double maxWidth = 1;
        if (featureList != null && featureList.size() > 0) {


            // Get max feature width
            for (BedPEFeature feature : featureList) {
                // Note -- don't cast these to an int until the range is checked.
                // could get an overflow.
                double pixelStart = ((feature.getStart() - origin) / locScale);
                double pixelEnd = ((feature.getEnd() - origin) / locScale);

                // If the any part of the feature fits in the
                // Track rectangle draw it
                if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {
                    if (pixelEnd - pixelStart > maxWidth) maxWidth = pixelEnd - pixelStart;
                }
            }
        }

        double f = trackRectangle.height / maxWidth;

        for (BedPEFeature feature : featureList) {

            // Note -- don't cast these to an int until the range is checked.
            // could get an overflow.
            double pixelStart = ((feature.getStart() - origin) / locScale);
            double pixelEnd = ((feature.getEnd() - origin) / locScale);

            // If the any part of the feature fits in the Track rectangle draw it
            if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {


                Color color = track.getColor();
                Graphics2D g = context.getGraphic2DForColor(color);

                int w = (int) (pixelEnd - pixelStart);
                if (w < 3) {
                    w = 3;
                    pixelStart--;
                }

                double h = Math.min(trackRectangle.height, f * w);
                double y = trackRectangle.y + trackRectangle.height - h;

                Arc2D.Double arcPath = new Arc2D.Double(pixelStart, y, w, 2 * h, 0, 180, Arc2D.OPEN);

                g.draw(arcPath);

            }
        }
    }


}
