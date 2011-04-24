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
public class PeakRenderer implements Renderer<Peak> {


    public void render(List<Peak> peakList, RenderContext context, Rectangle rect, Track track) {

        final double locScale = context.getScale();
        final double origin = context.getOrigin();

        int lastPX = -1;
        double pXMin = rect.getMinX();
        double pXMax = rect.getMaxX();

        float[] bgColorComps = context.getBackgroundColor().getColorComponents(null);
        float[] fgColorComps = track.getColor().getColorComponents(null);

        for (Peak peak : peakList) {

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
            if (PeakTrack.getShadeOption() == PeakTrack.ShadeOption.FOLD_CHANGE) {
                score = peak.getDynamicScore();
                if (score > 0 && score < 1) {
                    score = 1 / score;
                }
            }
            if (score < ((PeakTrack) track).getScoreThreshold()) continue;

            int top = rect.y;
            int h = rect.height;
            final float[] timeScores = peak.getTimeScores();
            if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {
                h = (rect.height) / (timeScores.length + 1);
            }
            drawScore(context, bgColorComps, fgColorComps, pX, dX, top, h, score);

            if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {

                for (int i = 0; i < timeScores.length; i++) {
                    top += h;
                    score = timeScores[i];
                    drawScore(context, bgColorComps, fgColorComps, pX, dX, top, h, score);

                }
            }

        }
    }

    private void drawScore(RenderContext context, float[] bgColorComps, float[] fgColorComps,
                           int pX, int dX, int top, int h, float score) {
        Color c = null;
        //if (peak.isDynamic()) {
        //    c = Color.red;
        //} else {
        float alpha = 1.0f;
        if (PeakTrack.getShadeOption() == PeakTrack.ShadeOption.SCORE) {
            // scale is 1 -> 100
            int shadeStep = (int) (score / 10);
            alpha = Math.max(0.1f, (Math.min(1.0f, shadeStep * 0.1f)));
        } else if (PeakTrack.getShadeOption() == PeakTrack.ShadeOption.FOLD_CHANGE) {
            // Scale is 3 -> 6
            if (score < 3) {
                alpha = 0.1f;
            } else {
                int shadeStep = (int) (score / .6);
                alpha = Math.max(0.1f, (Math.min(1.0f, shadeStep * 0.1f)));
            }
        }

        c = ColorUtilities.getCompositeColor(bgColorComps, fgColorComps, alpha);
        //}

        Graphics2D g = context.getGraphic2DForColor(c);

        g.fillRect(pX, top + 1, dX, h - 2);
    }
}
