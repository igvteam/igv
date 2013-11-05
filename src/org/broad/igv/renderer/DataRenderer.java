/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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

    /**
     * Render the provided scores. No border, scales, axes, or anything else
     * @param track
     * @param scores
     * @param context
     * @param arect
     */
    public abstract void renderScores(Track track, List<LocusScore> scores,
                                         RenderContext context, Rectangle arect);

}
