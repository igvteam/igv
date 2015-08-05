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
 * TrackRenderer.java
 *
 * Created on Sep 6, 2007, 10:07:39 AM
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

import org.broad.igv.PreferenceManager;
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
    public void renderScores(Track track, List<LocusScore> locusScores, RenderContext context, Rectangle arect) {


        Rectangle adjustedRect = calculateDrawingRect(arect);


        double origin = context.getOrigin();
        double locScale = context.getScale();

        Color posColor = track.getColor();
        Color negColor = track.getAltColor();

        Graphics2D gPos = context.getGraphic2DForColor(posColor);
        Graphics2D gNeg = context.getGraphic2DForColor(negColor);

        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
            gPos.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            gNeg.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        }


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
        gPos.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_DEFAULT);
        gNeg.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_DEFAULT);

    }
}
