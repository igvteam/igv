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
 * DataRenderer.java
 *
 * Created on November 27, 2007, 9:20 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */


package org.broad.igv.renderer;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 */
public abstract class DataRenderer implements Renderer<LocusScore> {

    private static Logger log = Logger.getLogger(DataRenderer.class);

    protected static final int AXIS_AREA_WIDTH = 60;
    protected static Color axisLineColor = new Color(255, 180, 180);

    /**
     * Render the track in the given rectangle.
     *
     * @param track
     * @param scores
     * @param context
     * @param rect
     */
    public void render(List<LocusScore> scores, RenderContext context, Rectangle rect, Track track) {

        if (scores != null) {
            // Prevent modification of the scores collection during rendering.  This collection
            // has caused concurrent modification exceptions.
            synchronized (scores) {
                renderScores(track, scores, context, rect);
                renderAxis(track, context, rect);
            }
        }
        renderBorder(track, context, rect);

    }

    /**
     * Render a border.  By default does nothing.
     *
     * @param track
     * @param context
     * @param rect
     */
    public void renderBorder(Track track, RenderContext context, Rectangle rect) {
    }

    /**
     * Render a Y axis.  By default does nothing.
     *
     * @param track
     * @param context
     * @param rect
     */
    public void renderAxis(Track track, RenderContext context, Rectangle rect) {
        PreferenceManager prefs = PreferenceManager.getInstance();

        // For now disable axes for all chromosome view
        if (context.getChr().equals(Globals.CHR_ALL)) {
            return;
        }
        if (prefs.getAsBoolean(PreferenceManager.CHART_DRAW_Y_AXIS))  {

            Rectangle axisRect = new Rectangle(rect.x, rect.y + 1, AXIS_AREA_WIDTH, rect.height);
            Graphics2D whiteGraphics = context.getGraphic2DForColor(Color.white);

            whiteGraphics.fillRect(axisRect.x, axisRect.y, axisRect.width, axisRect.height);

            Graphics2D axisGraphics = context.getGraphic2DForColor(axisLineColor);

            axisGraphics.drawLine(rect.x + AXIS_AREA_WIDTH, rect.y, rect.x + AXIS_AREA_WIDTH,
                    rect.y + rect.height);
        }


    }

    protected abstract void renderScores(Track track, List<LocusScore> scores,
                                         RenderContext context, Rectangle arect);

}
