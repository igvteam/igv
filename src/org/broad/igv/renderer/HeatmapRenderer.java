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

package org.broad.igv.renderer;

import org.broad.igv.PreferenceManager;
import org.broad.igv.data.rnai.RNAIGeneScore;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;

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
        boolean showAllFeatures = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_SHOW_ALL_HEATMAP);

        double maxX = rect.getMaxX();
        int minY = (int) rect.getMinY();
        int height = (int) rect.getHeight();
        int lastPEnd = 0;
        int lastPStart = 0;
        int lastW = 0;

        for (LocusScore score : scores) {
            if (lastPStart > maxX) {
                break;
            }

            // Note -- don't cast these to an int until the range is checked,
            // otherwise could get an overflow.
            float fStart = (float) ((score.getStart() - origin) / locScale);
            float fEnd = (float) ((score.getEnd() - origin) / locScale);
            // float fw = fEnd - fStart;
            int pStart = (int) fStart;
            int pEnd = (int) fEnd;

            int min;
            if(showAllFeatures) {
                min = 1;
            } else {
                min = Math.min(pStart - lastPEnd, 1);
            }

            int w = Math.max(min, pEnd - pStart);

            float dataY = track.logScaleData(score.getScore());
            Color graphColor = colorScale.getColor(dataY);

            if ((pStart + w) >= 0 && (lastPStart <= maxX)) {

                // Minimum width for DNA methylation
                if(track.getTrackType() == TrackType.DNA_METHYLATION) {
                    if(w < 6) {
                        int pMid = (pStart + pEnd) / 2;
                        pStart = pMid - 3;
                        w = 6;
                    }
                }

                // This test handles the rather pathological case where the previous feature was 1 pixel wide, and
                // the current feature overlaps it because of roundoff error when scaling.
                else if (pStart < lastPEnd && w > 1 && lastW == 1) {
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
                } // End special RNAi treagment
            }
            lastPStart = pStart;
            lastPEnd = pStart + w;
            lastW = w;
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
