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
import org.broad.igv.ui.FontManager;

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


    /**
     * Draw scale in top left of rectangle
     * @param range
     * @param context
     * @param arect
     */
    public static void drawScale(DataRange range, RenderContext context, Rectangle arect){
        if (range != null) {
            Graphics2D g = context.getGraphic2DForColor(Color.black);
            Font font = g.getFont();
            Font smallFont = FontManager.getFont(8);
            try {
                g.setFont(smallFont);
                String minString = range.getMinimum() == 0f ? "0" : String.format("%.3f", range.getMinimum());
                String fmtString = range.getMaximum() > 10 ? "%.0f" : "%.2f";
                String maxString = String.format(fmtString, range.getMaximum());
                String scale = "[" + minString + " - " + maxString + "]";
                g.drawString(scale, arect.x + 5, arect.y + 10);

            } finally {
                g.setFont(font);
            }
        }
    }

}
