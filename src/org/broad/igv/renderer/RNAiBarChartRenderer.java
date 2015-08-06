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
