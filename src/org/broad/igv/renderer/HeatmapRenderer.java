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

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.seg.Segment;
import org.broad.igv.data.rnai.RNAIGeneScore;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ColorUtilities;

import java.awt.*;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class HeatmapRenderer extends DataRenderer {




    public String getDisplayName() {
        return "Heatmap";
    }

    /**
     * Render the data track as a heat map.
     * <p/>
     * This method has gotten quite complicated,  most of it from the option to join adjacent
     * copy number segments.
     *
     * @param track
     * @param scores
     * @param context
     * @param rect
     */
    public void renderScores(Track track, List<LocusScore> scores, RenderContext context, Rectangle rect) {

        ContinuousColorScale colorScale = track.getColorScale();


        double origin = context.getOrigin();
        double locScale = context.getScale();

        Color bgColor = colorScale.getNoDataColor();
        context.getGraphic2DForColor(bgColor).fill(rect);


        double maxX = rect.getMaxX();
        int minY = (int) rect.getMinY();
        int height = (int) rect.getHeight();
        int lastPEnd = 0;
        int lastPStart = 0;
        int lastW = 0;
        Color lastColor = null;
        // These buffers are created here to be thread-safe
        float[] buffer1 = new float[3];
        float[] buffer2 = new float[3];

        for (LocusScore score : scores) {
            if (lastPStart > maxX) {
                break;
            }

            // Note -- don't cast these to an int until the range is checked,
            // otherwise could get an overflow.
            float fStart = (float) ((score.getStart() - origin) / locScale);
            float fEnd = (float) ((score.getEnd() - origin) / locScale);
            float fw = fEnd - fStart;
            int pStart = (int) fStart;
            int pEnd = (int) fEnd;
            int w = Math.max(1,  pEnd - pStart); 

            // if the width is < 1 pixel use alpha to mix color with last color, or background
            float dataY = track.logScaleData(score.getScore());
            Color graphColor = colorScale.getColor(dataY);
            if (fw < 1) {
                float alpha = Math.max(0.25f, fw);
                if(lastColor == null || pStart > lastPEnd) {
                    graphColor = ColorUtilities.getCompositeColor(graphColor.getColorComponents(buffer1), alpha);
                }
                else {
                    graphColor = ColorUtilities.getCompositeColor(lastColor.getColorComponents(buffer1),
                            graphColor.getColorComponents(buffer2), alpha);
                }
            } 

            if ((pStart + w) >= 0 && (lastPStart <= maxX)) {

                // This test handles the rather pathological case where the previous feature was 1 pixel wide, and
                // the current feature overlaps it because of roundoff error when scaling.
                if (pStart < lastPEnd && w > 1 && lastW == 1) {
                    pStart++;
                    w--;
                }

                // TODo The instanceof test is very very ugly.   Special RNAi treatment
                // Refactor  to generalize "confidence" for all datasets
                if (!Float.isNaN(dataY)) {
                    if (score instanceof RNAIGeneScore) {
                        RNAIGeneScore rnaiScore = (RNAIGeneScore) score;
                        if (rnaiScore.getConfidence() < 2) {
                            graphColor = getLowConfColor(context.getZoom());
                        }
                    }

                    Graphics2D g2D = context.getGraphic2DForColor(graphColor);
                    if (pStart < maxX) {
                        // Clip at edges
                        int pLeft = Math.max(rect.x, pStart);
                        int pRight = Math.min(rect.x + rect.width, pStart + w);
                        int adjustedW = pRight - pLeft;
                        g2D.fillRect(pLeft, minY, adjustedW, height);
                    }

                    // Segmented copy numbers (score type "segment") can be optionally joined
                    if (score instanceof Segment &&
                            PreferenceManager.getInstance().getAsBoolean(PreferenceManager.JOIN_ADJACENT_SEGMENTS_KEY) &&
                            lastColor != null && (pStart - lastPEnd) > 1) {

                        int midPoint = (pStart + lastPEnd) / 2;

                        context.getGraphic2DForColor(lastColor).fillRect(lastPEnd, minY,
                                midPoint - lastPEnd, height);
                        g2D.fillRect(midPoint, minY, pStart - midPoint, height);

                        // Cross hatch joined area  --- don't do for whole chromosome
                        if (!context.getChr().equals(Globals.CHR_ALL)) {
                            if (pStart - lastPEnd > 4 && height > 2) {
                                Color c = new Color(0.4f, 0.4f, 0.4f, 0.5f);
                                Graphics2D gLine = context.getGraphic2DForColor(c);
                                int midpoint = minY + height / 2;
                                if (height > 4) {
                                    gLine.drawLine(lastPEnd, minY + 3, lastPEnd, minY + height - 3);
                                    gLine.drawLine(pStart, minY + 3, pStart, minY + height - 3);
                                }
                                gLine.drawLine(lastPEnd, minY + height / 2, pStart - 1, minY + height / 2);
                            }
                        }
                    }
                } // End special RNAi treagment
            }
            lastPStart = pStart;
            lastPEnd = pStart + w;
            lastW = w;
            lastColor = graphColor;
        }
    }


    /**
     * Return the color indicating a low confdence score as a function of zoom level.  Currently
     * this is only used with RNAi data.
     *
     * @param zoom
     * @return
     */
    private static Color getLowConfColor(int zoom) {
        Color lowConfColor = lowConfColorCache.get(zoom);
        if (lowConfColor == null) {

            int value = 225 - Math.min(70, zoom * 10);
            lowConfColor = new Color(value, value, value);
            lowConfColorCache.put(zoom, lowConfColor);
        }
        return lowConfColor;
    }

    /**
     * A cache for "low confidence color", which is a function of zoom level.
     */
    static Map<Integer, Color> lowConfColorCache = new Hashtable();


}
