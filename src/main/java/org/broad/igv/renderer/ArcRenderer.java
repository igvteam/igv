package org.broad.igv.renderer;

import org.broad.igv.feature.IGVFeature;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.awt.geom.GeneralPath;
import java.util.List;

/**
 * Created by jrobinson on 9/17/15.
 */
public class ArcRenderer extends FeatureRenderer {

    @Override
    public void render(List<IGVFeature> featureList, RenderContext context, Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();
        Color ARC_COLOR = new Color(150, 50, 50);


        double maxWidth = 1;
        if (featureList != null && featureList.size() > 0) {

           // Get max feature width
            for (IGVFeature feature : featureList) {
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

        for (IGVFeature feature : featureList) {
            // Note -- don't cast these to an int until the range is checked.
            // could get an overflow.
            double pixelStart = ((feature.getStart() - origin) / locScale);
            double pixelEnd = ((feature.getEnd() - origin) / locScale);

            // If the any part of the feature fits in the
            // Track rectangle draw it
            if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {

                // Set color used to draw the feature

                Color color = feature.getColor();
                if(color == null) {
                    color = track.getColor();
                }
                Graphics2D g = context.getGraphic2DForColor(color);
                g.setFont(FontManager.getDefaultFont());


                int w = (int) (pixelEnd - pixelStart);
                if (w < 3) {
                    w = 3;
                    pixelStart--;
                }

                double h = Math.min(trackRectangle.height, f * w);
                double y = trackRectangle.y + trackRectangle.height - h;
                Arc2D.Double arcPath = new Arc2D.Double(pixelStart, y, w, 2*h, 0, 180, Arc2D.OPEN);


                g.draw(arcPath);

            }
        }
    }

}


