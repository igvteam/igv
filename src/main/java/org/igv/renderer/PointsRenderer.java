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
public class PointsRenderer extends XYPlotRenderer {

    int pointsSize = 2;

    public int getPointsSize() {
        return pointsSize;
    }

    public void setPointsSize(int pointsSize) {
        this.pointsSize = pointsSize;
    }

    protected void drawDataPoint(Color graphColor, int dx, int pX, int baseY, int pY, LocusScore score, RenderContext context) {
        if (context.getScale() < 4) {
            pointsSize = 4;
        } else {
            pointsSize = 2;
        }

        context.getGraphic2DForColor(graphColor).fillRect(pX, pY, dx, pointsSize);

    }


}
