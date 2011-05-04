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

package org.broad.igv.peaks;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.BarChartRenderer;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.util.ColorUtilities;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Apr 23, 2011
 */
public class PeakRenderer implements Renderer<LocusScore> {

    BarChartRenderer chartRenderer = new BarChartRenderer();
    public static final int SIGNAL_CHART_HEIGHT = 15;

    public void render(List<LocusScore> peakList, RenderContext context, Rectangle rect, Track track) {

        final double locScale = context.getScale();
        final double origin = context.getOrigin();

        int lastPX = -1;
        double pXMin = rect.getMinX();
        double pXMax = rect.getMaxX();

        float[] bgColorComps = context.getBackgroundColor().getColorComponents(null);
        float[] fgColorComps = track.getColor().getColorComponents(null);
        Graphics2D greyGraphics = context.getGraphic2DForColor(Color.lightGray);
        Graphics2D borderGraphics = context.getGraphic2DForColor(Color.black);

        int h = rect.height;
        int nTimePoints = ((PeakTrack) track).nTimePoints;
        if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {
            h = (rect.height - 2) / (nTimePoints + 1);
        }

        for (LocusScore ls : peakList) {

            Peak peak = (Peak) ls;

            int start = peak.getStart();
            int end = peak.getEnd();

            int pX = (int) ((start - origin) / locScale);
            int dX = (int) Math.max(2, (end - start) / locScale);

            if (pX + dX < pXMin) {
                continue;
            }
            if (pX > pXMax) {
                break;
            }

            float score = peak.getCombinedScore();
            if (PeakTrack.getColorOption() == PeakTrack.ColorOption.FOLD_CHANGE) {
                score = peak.getDynamicScore();
            }

            int top = rect.y + 2;

            if (PeakTrack.getRenderOption() == PeakTrack.RenderOption.FEATURE) {
                drawScore(context, bgColorComps, fgColorComps, pX, dX, top, h, score, PeakTrack.getColorOption());
            }
            if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {
                float[] timeScores = peak.getTimeScores();
                for (int i = 0; i < timeScores.length; i++) {
                    top += h;
                    score = timeScores[i];
                    drawScore(context, bgColorComps, fgColorComps, pX, dX, top, h, score, PeakTrack.ColorOption.SCORE);
                }

            }

        }

        if (PeakTrack.getRenderOption() == PeakTrack.RenderOption.SIGNAL) {
            String chr = context.getChr();
            int start = (int) context.getOrigin();
            int end = (int) context.getEndLocation();
            int zoom = context.getZoom();

            List<LocusScore> signals =
                    ((PeakTrack) track).signalSource.getSummaryScoresForRange(chr, start, end, zoom);
            Rectangle signalRect = new Rectangle(rect.x, rect.y + 1, rect.width, h - 1);

            chartRenderer.render(signals, context, signalRect, track);
        }

        if (PeakTrack.getRenderOption() == PeakTrack.RenderOption.CHART) {
            Rectangle signalRect = new Rectangle(rect.x, rect.y + 1, rect.width, h - 1);
            chartRenderer.render(peakList, context, signalRect, track);
        }

        if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {
            borderGraphics.drawLine(rect.x, rect.y, rect.x + rect.width, rect.y);
            borderGraphics.drawLine(rect.x, rect.y + rect.height, rect.x + rect.width, rect.y + rect.height);
            greyGraphics.drawLine(rect.x, rect.y + h + 1, rect.x + rect.width, rect.y + h + 1);


        }
    }


    static float[] blueComps = Color.blue.getColorComponents(null);
    static float[] redComps = Color.red.getColorComponents(null);

    private void drawScore(RenderContext context, float[] bgColorComps, float[] fgColorComps,
                           int pX, int dX, int top, int h, float score, PeakTrack.ColorOption option) {
        Color c = null;
        //if (peak.isDynamic()) {
        //    c = Color.red;
        //} else {
        float alpha = 1.0f;
        if (option == PeakTrack.ColorOption.SCORE) {
            // scale is 1 -> 100
            int shadeStep = (int) (score / 10);
            alpha = Math.max(0.1f, (Math.min(1.0f, shadeStep * 0.1f)));

        } else if (option == PeakTrack.ColorOption.FOLD_CHANGE) {

            // Scale is 1.58 -> 3.3 log2
            if (Math.abs(score) < 1.58) {
                alpha = 1f;
                fgColorComps = Color.gray.getColorComponents(null);
            } else {
                // vary alpha from .3 -> 1 over range 3.3 -> 1.58 in steps of 10
                int shadeStep = (int) ((Math.abs(score) - 1.58) / (3.3 - 1.58) * 10);
                alpha = Math.max(0.3f, (Math.min(1f, 0.3f + shadeStep * 0.1f)));
                fgColorComps = (score < 0 ? blueComps : redComps);
            }
        }

        c = ColorUtilities.getCompositeColor(bgColorComps, fgColorComps, alpha);
        //}

        Graphics2D g = context.getGraphic2DForColor(c);

        g.fillRect(pX, top + 1, dX, h - 2);
    }
}
