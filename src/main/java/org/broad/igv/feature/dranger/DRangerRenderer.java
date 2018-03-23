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
