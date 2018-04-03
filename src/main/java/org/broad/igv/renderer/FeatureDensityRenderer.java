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

import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class FeatureDensityRenderer extends DataRenderer {

    public String getDisplayName() {
        return "Feature Density";
    }

    /**
     * Render the data track as a bar chart.
     */
    public void renderScores(Track track, List<LocusScore> scores, RenderContext context, Rectangle rect) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        double maxValue = track.getDataRange().getMaximum();

        double yScaleFactor = rect.getHeight() / maxValue;

        int pX = -1;

        Graphics2D g = context.getGraphics();

        for (LocusScore score : scores) {
            pX = (int) ((score.getStart() - origin) / locScale);

            if (pX >= 0) {

                // Plot as density in counts per MB.  Assuming bin width is 1 pixel,
                // bin width in BP is locScale * 1.  So density in bp is score / locScale.
                //    TODO -- use actual bin width
                float dataY = (float) ((score.getScore() * 1000000) / locScale);

                if (Float.isNaN(dataY)) {
                    g.setColor(Color.LIGHT_GRAY.brighter());
                    g.drawLine(pX, (int) rect.getY(), pX, (int) rect.getMaxY());

                } else {
                    g.setColor(Color.CYAN.darker());
                    double scaledY = dataY * yScaleFactor;
                    int pY = (int) Math.max(rect.getY(), (rect.getMaxY() - scaledY));
                    g.drawLine(pX, (int) rect.getMaxY(), pX, pY);
                }
                if (pX > rect.getMaxX()) {
                    break;
                }


                // Draw a single dividing line along the bottom or the rect
                g.setColor(Color.BLACK);
                g.drawLine((int) rect.getX(), (int) rect.getMaxY(), (int) rect.getMaxX(), (int) rect.getMaxY());
            }
        }

    }

}


