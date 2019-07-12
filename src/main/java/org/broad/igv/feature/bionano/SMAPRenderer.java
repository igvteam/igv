/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.feature.bionano;

import org.broad.igv.feature.IGVFeature;
import org.broad.igv.renderer.FeatureRenderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.util.List;


/**
 * Created by jrobinson on 9/17/15.
 */
public class SMAPRenderer extends FeatureRenderer {

    static final int BLOCK_HEIGHT = 14;

    @Override
    public void render(List<IGVFeature> featureList, RenderContext context, Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();


        for (IGVFeature feature : featureList) {
            // Note -- don't cast these to an int until the range is checked.
            // could get an overflow.
            double pixelStart = Math.max(trackRectangle.x, ((feature.getStart() - origin) / locScale));
            double pixelEnd = Math.min(trackRectangle.x + trackRectangle.width, ((feature.getEnd() - origin) / locScale));

            if (pixelEnd < trackRectangle.x) continue;
            if (pixelStart > (trackRectangle.x + trackRectangle.width)) break;

            Color color = feature.getColor();
            if (color == null) {
                color = track.getColor();
            }
            Graphics2D g = context.getGraphic2DForColor(color);
            int y = trackRectangle.y + trackRectangle.height / 2;

            if (feature instanceof SMAPPairedFeature) {

                SMAPFeature feature1 = ((SMAPPairedFeature) feature).feature1;
                SMAPFeature feature2 = ((SMAPPairedFeature) feature).feature2;
                int x1 =  (int) Math.max(trackRectangle.x, ((feature1.getEnd() - origin) / locScale));
                int x2 =  (int) Math.min(trackRectangle.x + trackRectangle.width, ((feature2.getStart() - origin) / locScale));

                g.drawLine(x1, y, x2, y);
                drawSMapFeature(g, trackRectangle, origin, locScale, y, feature1);
                drawSMapFeature(g, trackRectangle, origin, locScale, y, feature2);

            } else {
                drawSMapFeature(g, trackRectangle, origin, locScale, y, feature);
            }
        }
    }

    public void drawSMapFeature(Graphics2D g, Rectangle trackRectangle, double origin, double locScale, int y, IGVFeature feature) {

        double pixelStart = Math.max(trackRectangle.x, ((feature.getStart() - origin) / locScale));
        double pixelEnd = Math.min(trackRectangle.x + trackRectangle.width, ((feature.getEnd() - origin) / locScale));

        // If the any part of the feature fits in the
        // Track rectangle draw it
        if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {

            // Set color used to draw the feature


            int w = (int) (pixelEnd - pixelStart);
            if (w < 3) {
                w = 3;
                pixelStart--;
            }
            g.fillRect((int) pixelStart, y - BLOCK_HEIGHT / 2, w, BLOCK_HEIGHT);

        }
    }

}


