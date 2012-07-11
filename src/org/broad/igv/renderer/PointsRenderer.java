/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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

import org.broad.igv.track.RenderContext;

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

    protected void drawDataPoint(Color graphColor, int dx, int pX, int baseY, int pY, RenderContext context) {
        if (context.getScale() < 4) {
            pointsSize = 4;
        } else {
            pointsSize = 2;
        }

        context.getGraphic2DForColor(graphColor).fillRect(pX, pY, dx, pointsSize);

    }


}
