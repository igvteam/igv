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

import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.GisticScore;
import org.broad.igv.track.GisticTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class GisticTrackRenderer {

    public void setOverlayMode(boolean mode) {
        // Ignore for now
    }

    public void render(Track track, RenderContext context, Rectangle rect) {

        // I only know how to render gistic tracks
        if (!(track instanceof GisticTrack)) {
            return;
        }
        plotScores(track, context, rect);
    }

    protected void plotScoresOn(List<GisticScore> scores, Graphics2D g2D, Rectangle rect, int xEnd,
                                double scale, RenderContext vc, int xStart, int yStart) {

        double origin = vc.getOrigin();
        double locScale = vc.getScale();

        int lastY = yStart;
        double lastX = xStart;
        for (GisticScore score : scores) {
            double scaledY = scale * transform(score.getQValue());
            int y = (int) (rect.getMaxY() - scaledY);

            // Compute X values in double precision to prevent overflow at extreme
            // scales.  Cast to ints after range is checked.
            double x1 = ((score.getStart() - 1 - origin) / locScale);
            double x2 = ((score.getEnd() - origin) / locScale);
            if (x2 > rect.getX() && lastX < rect.getMaxX()) {
                g2D.drawLine((int) lastX, lastY, (int) x1, y);
                g2D.drawLine((int) x1, y, (int) x2, y);
            }
            lastY = y;
            lastX = x2;
        }
        if (xEnd > rect.getX() && lastX < rect.getMaxX()) {
            g2D.drawLine((int) lastX, lastY, xEnd, yStart);
        }
    }

    // TODO -- need to combine the amp and del regions with a parameter.  Remove
    // duplicated code.
    private void plotScores(Track track, RenderContext context, Rectangle rect) {

        GisticTrack gisticTrack = (GisticTrack) track;
        String chr = context.getChr();

        float maxQValue = transform(track.getDataRange().getMaximum());
        if (maxQValue > 0) {
            double scale = (rect.getHeight() - 5) / maxQValue;

            // TODO use offset
            int xStart = (int) rect.getX();
            int yStart = (int) rect.getMaxY();

            //long chromosomeEnd = offset + Genome.getInstance().getChromosomeLength(chr);
            int xEnd = (int) rect.getWidth(); //vc.getPy(chromosomeEnd) - 1;

            List<GisticScore> scores = gisticTrack.getAmpScores(chr);
            if (scores != null) {
                Graphics2D g2D = context.getGraphic2DForColor(Color.RED);
                plotScoresOn(scores, g2D, rect, xEnd, scale, context, xStart, yStart);
            }

            scores = gisticTrack.getDelScores(chr);
            if (scores != null) {
                Graphics2D g2D = context.getGraphic2DForColor(Color.BLUE);
                plotScoresOn(scores, g2D, rect, xEnd, scale, context, xStart, yStart);
            }
        }
    }

    public float transform(double value) {
        return (float) Math.log(1 + value);
    }

    /**
     * Render a border.  By default does nothing.
     */
    public void renderBorder(Track track, RenderContext context, Rectangle rect) {
        context.getGraphic2DForColor(Color.gray).drawRect(rect.x, rect.y, rect.width, rect.height);
    }

    /**
     * Render a Y axis.  TODO -- implementation
     *
     * @param track
     * @param context
     * @param rect
     */
    public void renderAxis(Track track, RenderContext context, Rectangle rect) {

    }
}
