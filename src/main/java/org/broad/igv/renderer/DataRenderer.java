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

import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;

import java.awt.*;
import java.util.List;

import static org.broad.igv.prefs.Constants.CHART_DRAW_Y_AXIS;

/**
 * @author jrobinso
 */
public abstract class DataRenderer implements Renderer<LocusScore> {

    private static Logger log = LogManager.getLogger(DataRenderer.class);

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
        IGVPreferences prefs = PreferencesManager.getPreferences();

        // For now disable axes for all chromosome view
        if (context.getChr().equals(Globals.CHR_ALL)) {
            return;
        }
        if (prefs.getAsBoolean(CHART_DRAW_Y_AXIS)) {

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
     *
     * @param track
     * @param scores
     * @param context
     * @param arect
     */
    public abstract void renderScores(Track track, List<LocusScore> scores,
                                      RenderContext context, Rectangle arect);


    /**
     * Draw scale in top left of rectangle
     *
     * @param range
     * @param context
     * @param arect
     */
    public static void drawScale(DataRange range, RenderContext context, Rectangle arect) {
        if (range != null && context.multiframe == false) {
            Graphics2D g = (Graphics2D) context.getGraphics().create();
            if(Globals.isDarkMode()) {
                g.setColor(Color.WHITE);
            } else {
                g.setColor(Color.BLACK);
            }
            Font smallFont = FontManager.getFont(8);
            try {
                g.setFont(smallFont);
                String minString = range.getMinimum() == 0f ? "0" : String.format("%.3f", range.getMinimum());
                String fmtString = range.getMaximum() > 10 ? "%.0f" : "%.2f";
                String maxString = String.format(fmtString, range.getMaximum());
                String scale = "[" + minString + " - " + maxString + "]";
                g.drawString(scale, arect.x + 5, arect.y + 10);

            } finally {
                g.dispose();
            }
        }
    }

}
