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


package org.broad.igv.renderer;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.data.rnai.RNAIGeneScore;
import org.broad.igv.data.rnai.RNAIHairpinValue;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class RNAiBarChartRenderer extends XYPlotRenderer {


    /**
     * Render the track in the given rectangle.
     *
     * @param track
     * @param locusScores
     * @param context
     * @param trackRect
     */
    @Override
    public void renderScores(Track track, List<LocusScore> locusScores, RenderContext context,
                             Rectangle trackRect) {


        // Experimental -- adjust the rectangle to maintain a buffer (space) between tracks.  The
        // settings (20% with 10 pixel maximum) should be adjustable.
        double buffer = Math.min(trackRect.getHeight() * 0.2, 10);
        Rectangle drawingRect = new Rectangle(trackRect);
        drawingRect.y = (int) (trackRect.getY() + buffer);
        drawingRect.height = (int) (trackRect.height - (drawingRect.y - trackRect.getY()));

        // The starting location (origin) in bp (base pairs) of the current view
        double origin = context.getOrigin();

        // The locaction scale  in bp / pixel of the current view
        double locScale = context.getScale();

        // The axis definition defines the scale of the y axis (minimum, base, and maximum values)
        DataRange axisDefinition = track.getDataRange();

        // Calculate the certical (Y) position in pixels of the base value.
        int baseY = computeYPixelValue(drawingRect, axisDefinition, axisDefinition.getBaseline());

        for (LocusScore score : locusScores) {

            // We'll have to case this score to an RNAiGeneScore here for now.  This is not
            // particularly safe and should be redesigned to eliminate the need to cache.
            RNAIGeneScore rnaiScore = (RNAIGeneScore) score;

            // Note -- don't cast these to an int until the range is checked.
            // could get an overflow.
            double pX = ((score.getStart() - origin) / locScale);
            double dx = Math.ceil((score.getEnd() - score.getStart()) / locScale) + 1;

            // Don't bother drawing if the data is not in range
            if ((pX + dx) >= 0 && (pX <= drawingRect.getMaxX())) {
                float dataY = score.getScore();

                // Draw the gene score
                if (!Float.isNaN(dataY)) {
                    int pY = computeYPixelValue(drawingRect, axisDefinition, dataY);
                    drawDataPoint(new Color(57, 67, 201), (int) dx, (int) pX, baseY, pY,
                            context);
                }

                // Draw the hairpin values, if any
                if (rnaiScore.getHairpinValues() != null) {
                    // my blue color new Color(195, 211, 237)
                    for (RNAIHairpinValue hValue : rnaiScore.getHairpinValues()) {
                        int scoreY = computeYPixelValue(drawingRect,
                                axisDefinition,
                                hValue.getScoreMean());
                        context.getGraphic2DForColor(new Color(95, 120, 200).brighter()).drawLine(
                                (int) pX, scoreY, ((int) pX + (int) dx) - 1, scoreY);
                    }
                }
            }
        }
    }

    @Override

    /**
     * Render a data point,  in this case as a bar.
     */
    protected void drawDataPoint(Color graphColor, int dx, int pX, int baseY, int pY,
                                 RenderContext context) {

        if (pY > baseY) {
            context.getGraphic2DForColor(graphColor).fillRect(pX, baseY, dx, pY - baseY);
        } else {
            context.getGraphic2DForColor(graphColor).fillRect(pX, pY, dx, baseY - pY);
        }

        // }
    }


    public String getDisplayName() {
        return "RNAi BarChart";
    }
}
