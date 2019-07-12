package org.broad.igv.feature.dsi;

import htsjdk.tribble.Feature;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.util.List;

/**
 * Created by jrobinson on 7/19/16.
 */
public class DSIRenderer implements Renderer<Feature> {


    static final Color uColor = new Color(26, 204, 255);
    static final Color mColor = new Color(153, 253, 153);
    static final Color pColor = new Color(50, 200, 50);
    static final Color fColor = new Color(253, 102, 8);
    static final Color naColor = Color.lightGray;


    double max = 20.0;

    @Override
    public void render(List<Feature> featureList, RenderContext context, Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        if (featureList != null && featureList.size() > 0) {

            for (Feature f : featureList) {

                if (!(f instanceof DSIFeature)) continue;   // Don't know how to draw anything else

                DSIFeature feature = (DSIFeature) f;

                // Note -- don't cast these to an int until the range is checked.
                // could get an overflow.
                double pixelStart = ((feature.getStart() - origin) / locScale);
                double pixelEnd = ((feature.getEnd() - origin) / locScale);

                // If the any part of the feature fits in the
                // Track rectangle draw it
                if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {

                    int w = (int) (pixelEnd - pixelStart);
                    if (w < 1) {
                        w = 1;
                    }
                    else if(w > 5) {
                        w -= 2;
                        pixelStart += 1;
                    }


                    int totalHeight = trackRectangle.height; //((feature.total / max) * trackRectangle.getHeight());
                    int bottom = trackRectangle.y + trackRectangle.height;
                    int top = trackRectangle.y;
                    int x = (int) pixelStart;

                        Color color = naColor;
                        Graphics2D g2D = context.getGraphic2DForColor(color);
                        g2D.fillRect(x, top, w, totalHeight);

                    if (feature.u != Integer.MIN_VALUE) {

                        double total = feature.total;

                        g2D = context.getGraphic2DForColor(uColor);
                        int uHeight = (int) Math.round((feature.u / total) * totalHeight);
                        g2D.fillRect(x, top, w, uHeight);
                        top += uHeight;

                        g2D = context.getGraphic2DForColor(mColor);
                        int mHeight = (int) Math.round((feature.m / total) * totalHeight);
                        g2D.fillRect(x, top, w, mHeight);
                        top += mHeight;

                        g2D = context.getGraphic2DForColor(pColor);
                        int pHeight = (int) Math.round((feature.p / total) * totalHeight);
                        g2D.fillRect(x, top, w, pHeight);
                        top += pHeight;

                        g2D = context.getGraphic2DForColor(fColor);
                        int fHeight = totalHeight - top;
                        g2D.fillRect(x, top, w, fHeight);
                    }


                }
            }

        }
    }
}
