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


