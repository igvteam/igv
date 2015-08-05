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
