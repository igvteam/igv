/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import org.broad.igv.PreferenceManager;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.RenderContext;
import org.broad.tribble.Feature;

import java.awt.*;
import java.util.List;

/**
 * Just draws a box.
 */
public class BasicTribbleRenderer extends FeatureRenderer<Feature> {
    @Override
    public void renderFeatures(List<Feature> featureList, RenderContext context,
                               Rectangle trackRectangle, FeatureTrack track) {

        double origin = context.getOrigin();
        double scale = context.getScale();

        if (featureList != null && featureList.size() > 0) {

            Graphics2D g2d = context.getGraphic2DForColor(track.getColor());
            for (Feature feature : featureList) {

                double pixelStart = ((feature.getStart() - 1 - origin) / scale);
                double pixelEnd = ((feature.getEnd() - origin) / scale);

                // If the any part of the feature fits in the
                // Track rectangle draw it
                if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {

                    int boxWidth = (int) Math.max(1, (pixelEnd - pixelStart));
                    int boxHeight = (int) Math.max(1, trackRectangle.getHeight() - 2);
                    int boxY = (int) (trackRectangle.getY() + (trackRectangle.getHeight() - boxHeight) / 2);


                    Rectangle boxRect = new Rectangle((int) pixelStart, boxY, boxWidth, boxHeight);

                     g2d.fill(boxRect);
                }

                // Features are sorted by position, so if the start is beyond the track rectangle we're done
                if(pixelStart > trackRectangle.getMaxX()) {
                    break;
                }
            }
        }
    }
}
