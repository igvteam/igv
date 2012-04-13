/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
