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
 * TrackRenderer.java
 *
 * Created on Sep 6, 2007, 10:07:39 AM
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class LineplotRenderer extends XYPlotRenderer {

    /**
     * Render the track in the given rectangle.
     *
     * @param track
     * @param locusScores
     * @param context
     * @param arect
     */
    @Override
    public void renderScores(Track track, List<LocusScore> locusScores, RenderContext context,
                             Rectangle arect) {


        Rectangle adjustedRect = calculateDrawingRect(arect);


        double origin = context.getOrigin();
        double locScale = context.getScale();

        Color posColor = track.getColor();
        Color negColor = track.getAltColor();

        Graphics2D gPos = context.getGraphic2DForColor(posColor);
        gPos.setRenderingHint(RenderingHints.KEY_ANTIALIASING, // Anti-alias!
                RenderingHints.VALUE_ANTIALIAS_ON);
        Graphics2D gNeg = context.getGraphic2DForColor(negColor);
        gNeg.setRenderingHint(RenderingHints.KEY_ANTIALIASING, // Anti-alias!
                RenderingHints.VALUE_ANTIALIAS_ON);


        // Get the Y axis definition, consisting of minimum, maximum, and base value.  Often
        // the base value is == min value which is == 0.
        DataRange axisDefinition = track.getDataRange();
        float maxValue = (float) axisDefinition.getMaximum();
        float baseValue = (float) axisDefinition.getBaseline();
        float minValue = (float) axisDefinition.getMinimum();

        // Calculate the Y scale factor.
        double yScaleFactor = adjustedRect.getHeight() / (maxValue - minValue);

        int lastPx = 0;
        int lastPy = Integer.MIN_VALUE;
        for (LocusScore score : locusScores) {
            // Note -- don't cast these to an int until the range is checked.
            // could get an overflow.
            double x = ((score.getStart() - origin) / locScale);
            double dx = (score.getEnd() - score.getStart()) / locScale;

            float dataY = score.getScore();

            // Compute the pixel y location.  
            double y = adjustedRect.getY() + (maxValue - dataY) * yScaleFactor;

            if ((x + dx < 0 || lastPy == Integer.MIN_VALUE)) {
                // Offscreen.  Just record the points
                lastPx = (int) (x + dx);
                lastPy = (int) (y);
                continue;
            }


            // Missing data in a dataset is signifed by NaN.  Just skip these.
            if (!Float.isNaN(dataY)) {
                int pX = (int) x;
                int pY = (int) y;
                double slope = ((double) pY - lastPy) / (pX - lastPx);

                int clippedLastPX = lastPx;
                int clippedLastPY = lastPy;
                if (lastPy < adjustedRect.y || lastPy > (adjustedRect.y + adjustedRect.height)) {
                    clippedLastPY = (lastPy < adjustedRect.y) ? adjustedRect.y : adjustedRect.y + adjustedRect.height;
                    clippedLastPX = lastPx + (int) ((clippedLastPY - lastPy) / slope);
                }

                int clippedPX = pX;
                int clippedPY = pY;

                // Clip and interpolate if neccessary
                if (pY < adjustedRect.y || pY > (adjustedRect.y + adjustedRect.height)) {
                    clippedPY = (pY < adjustedRect.getMinY()) ? adjustedRect.y : adjustedRect.y + adjustedRect.height;
                    clippedPX = lastPx + (int) ((clippedPY - lastPy) / slope);
                }

                Graphics2D g = (dataY >= 0) ? gPos : gNeg;

                g.drawLine(clippedLastPX, clippedLastPY, clippedPX, clippedPY);


                if (dx >= 1 && (clippedPY == pY)) {
                    g.drawLine(pX, clippedPY, (int) (pX + dx), clippedPY);
                }

                lastPx = (int) (pX + dx);
                lastPy = (int) pY;

                if (lastPx > adjustedRect.getMaxX()) {
                    break;
                }


            }
        }

    }
}
