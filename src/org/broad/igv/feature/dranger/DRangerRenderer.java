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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.dranger;

import org.broad.igv.feature.IGVFeature;
import org.broad.igv.renderer.FeatureRenderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class DRangerRenderer extends FeatureRenderer {

    //TODO -- this is a copy from MutationRenderer.  Refactor to move this up
    // the class hierarchy.

    public void render(List<IGVFeature> featureList, RenderContext context,Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        if (featureList != null && featureList.size() > 0) {
            
            for (IGVFeature feature : featureList) {

                double pixelStart = ((feature.getStart() - origin) / locScale);
                double pixelEnd = ((feature.getEnd() - origin) / locScale);

                if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {
                    renderFeature(feature, context, trackRectangle);
                }
            }

        }
    }

    private void renderFeature(IGVFeature feature, RenderContext context, Rectangle rect) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        DRangerFeature dr = (DRangerFeature) feature;
        int x = (int) ((dr.getStart() - origin) / locScale);
        int length = dr.getEnd() - dr.getStart();
        int h = (int) Math.min(12, rect.getHeight());
        int w = (int) Math.ceil(length / locScale);
        int y = (int) Math.max(rect.getY(), (rect.getY() + (rect.getHeight() - h) / 2));

        Graphics2D g = context.getGraphic2DForColor(feature.getColor());
        g.fillRect(x, y, w, h);


    }


}
