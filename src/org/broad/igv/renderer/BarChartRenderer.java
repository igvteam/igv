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

import org.broad.igv.track.RenderContext;

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
    protected void drawDataPoint(Color graphColor, int dx, int pX, int baseY, int pY, RenderContext context) {
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
