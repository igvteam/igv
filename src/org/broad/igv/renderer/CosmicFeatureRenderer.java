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

package org.broad.igv.renderer;

import org.apache.log4j.Logger;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.color.ColorTable;
import org.broad.igv.ui.FontManager;
import org.broad.igv.PreferenceManager;

import java.awt.*;

public class CosmicFeatureRenderer extends FeatureRenderer {

    private static Logger log = Logger.getLogger(CosmicFeatureRenderer.class);
    ColorTable colorScheme;

    public CosmicFeatureRenderer() {
        colorScheme = PreferenceManager.getInstance().getMutationColorScheme();
    }

    public String getDisplayName() {
        return "Mutation";
    }


    /**
     * Note:  assumption is that featureList is sorted by start position.
     */
    public void render(java.util.List<IGVFeature> featureList, RenderContext context,
                       Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        if (featureList != null && featureList.size() > 0) {


            Rectangle lastRect = null;
            for (IGVFeature feature : featureList) {

                // Note -- don't cast these to an int until the range is checked.
                // could get an overflow.
                double pixelStart = Math.round((feature.getStart() - origin) / locScale);
                double pixelEnd = Math.round((feature.getEnd() - origin) / locScale);

                // If the any part of the feature fits in the
                // Track rectangle draw it
                if (pixelEnd >= trackRectangle.getX() &&
                        pixelStart <= trackRectangle.getMaxX()) {

                    // Set color used to draw the feature


                    int width = (int) pixelEnd - (int) pixelStart;
                    if (width < 3) {
                        width = 3;
                    }

                    int mutHeight = (int) Math.max(1, trackRectangle.getHeight() - 2);
                    int mutY = (int) (trackRectangle.getY() + (trackRectangle.getHeight() - mutHeight) / 2);

                    Rectangle mutRect = new Rectangle((int) pixelStart, mutY, width, mutHeight);

                    Color color = colorScheme.get(((IGVFeature) feature).getType());


                    Graphics2D g = context.getGraphic2DForColor(color);
                    Font font = FontManager.getDefaultFont();
                    g.setFont(font);

                    g.fill(mutRect);
                    if (lastRect != null && mutRect.intersects(lastRect)) {
                        // Indicate overlapping mutations
                        Graphics2D g2 = context.getGraphic2DForColor(Color.BLACK);
                        g2.draw(mutRect);

                    }


                    lastRect = mutRect;
                }
            }
        }
    }

}
