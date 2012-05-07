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
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.FontManager;

import java.awt.*;
import java.util.List;

public class MutationRenderer extends FeatureRenderer {

    private static Logger log = Logger.getLogger(MutationRenderer.class);

    public String getDisplayName() {
        return "Mutation";
    }

    /**
     * Note:  assumption is that featureList is sorted by start position.
     */
    public void render(List<IGVFeature> featureList, RenderContext context, Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        if (featureList != null && featureList.size() > 0) {


            Rectangle lastRect = null;

            boolean colorOverlay = IGV.getInstance().getSession().getColorOverlay();
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
                    Graphics2D g = context.getGraphic2DForColor(color);
                    g.setFont(FontManager.getDefaultFont());


                    int w = (int) (pixelEnd - pixelStart);
                    if (w < 3) {
                        w = 3;
                        pixelStart--;
                    }

                    int mutHeight = (int) Math.max(1, trackRectangle.getHeight() - 2);
                    int mutY = (int) (trackRectangle.getY() + (trackRectangle.getHeight() - mutHeight) / 2);


                    Rectangle mutRect = new Rectangle((int) pixelStart, mutY, w, mutHeight);

                    Graphics2D gRect = (colorOverlay ? g : context.getGraphic2DForColor(Color.BLACK));
                    gRect.draw(mutRect);
                    mutRect.x--;
                    mutRect.width += 2;
                    gRect.draw(mutRect);


                    lastRect = mutRect;
                }
            }

        }
    }
}
