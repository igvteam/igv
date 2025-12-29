/*
 * TrackRenderer.java
 *
 * Created on Sep 6, 2007, 10:07:39 AM
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.renderer;

import org.igv.feature.LocusScore;
import org.igv.track.RenderContext;

import java.awt.*;

/**
 * @author jrobinso
 */
public class BarChartRenderer extends XYPlotRenderer {

    @Override
    public String getDisplayName() {
        return "Bar Chart";
    }

    /**
     * Render the data track as a bar chart.
     */
    @Override
    protected void drawDataPoint(Color graphColor, int dx, int pX, int baseY, int pY, LocusScore score, RenderContext context) {
        //if (pY != baseY) {
        if (dx <= 1) {
            context.getGraphic2DForColor(graphColor).drawLine(pX, baseY, pX, pY);
        } else {
            //Graphics2D outlineContext = context.getGraphic2DForColor(Color.lightGray);
            if (pY > baseY) {
                context.getGraphic2DForColor(graphColor).fillRect(pX, baseY, dx, pY - baseY);
                //    outlineContext.drawRect(pX, baseY, dx, pY - baseY);

            } else {
                context.getGraphic2DForColor(graphColor).fillRect(pX, pY, dx, baseY - pY);
                //    outlineContext.drawRect(pX, pY, dx, baseY - pY);
            }
        }
        //}
    }
}
