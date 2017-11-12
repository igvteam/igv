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
* LocationPanel.java
*
* Created on September 11, 2007, 2:29 PM
*
* To change this template, choose Tools | Template Manager
* and open the template in the editor.
*/
package org.broad.igv.ui.javafx.panel;

//~--- non-JDK imports --------------------------------------------------------

import com.sun.javafx.tk.FontMetrics;
import com.sun.javafx.tk.Toolkit;
import javafx.geometry.Rectangle2D;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.ui.javafx.ResizableCanvas;
import org.broad.igv.ui.javafx.renderer.CytobandRenderer;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.List;

public class CytobandPane extends ResizableCanvas {

    private static double fontHeight = 10.0;
    private static double bandHeight = 10.0;
    private static String fontFamilyName = "Lucida Sans";
    private boolean isDragging = false;

    private double viewOrigin;
    private double viewEnd;

    private double cytobandScale;
    private ReferenceFrame frame;
    private Rectangle2D currentRegionRect;
    private CytobandRenderer cytobandRenderer;
    private List<Cytoband> currentCytobands;

    public CytobandPane(ReferenceFrame frame) {
        this(frame, true);
    }


    public CytobandPane(ReferenceFrame frame, boolean mouseable) {

        this.frame = frame;
        viewOrigin = frame.getOrigin();
        viewEnd = frame.getEnd();

        Font font = Font.font(fontFamilyName, FontWeight.BOLD, fontHeight);
        getCanvas().getGraphicsContext2D().setFont(font);

        // NOTE: use of internal com.sun.javafx APIs!
        FontMetrics fontMetrics = Toolkit.getToolkit().getFontLoader().getFontMetrics(font);

        cytobandRenderer = new CytobandRenderer(fontMetrics);
    }

    @Override
    protected void layoutChildren() {
        super.layoutChildren();

        if (frame.getChrName().equals(Globals.CHR_ALL) || getWidth() < 10) {
            return;
        }

        Chromosome chromosome = getReferenceFrame().getChromosome();
        if (chromosome == null) {
            return;
        }

        currentCytobands = chromosome.getCytobands();
        if (currentCytobands == null) {
            return;
        }

        GraphicsContext graphicContext = getCanvas().getGraphicsContext2D();

        int dataPanelWidth = frame.getWidthInPixels();
        Rectangle2D cytoRect = new Rectangle2D(0, 10, dataPanelWidth, bandHeight);

        cytobandRenderer.draw(currentCytobands, graphicContext, cytoRect, frame);

        int chromosomeLength = getReferenceFrame().getMaxCoordinate();
        cytobandScale = ((double) chromosomeLength) / dataPanelWidth;

        // The test is true if we are zoomed in
        if (getReferenceFrame().getZoom() > 0) {

            double origin = isDragging ? viewOrigin : getReferenceFrame().getOrigin();
            double end = isDragging ? viewEnd : getReferenceFrame().getEnd();

            double pixelStart = origin / cytobandScale;
            double pixelEnd = end / cytobandScale;
            double pixelSpan = Math.max(0, pixelEnd - pixelStart);

            // Draw Cytoband current region viewer
            double height = cytoRect.getHeight();
            graphicContext.setStroke(Color.RED);

            double y = cytoRect.getMinY() + CytobandRenderer.CYTOBAND_Y_OFFSET;
            currentRegionRect = new Rectangle2D(pixelStart - 2, y, pixelSpan + 4, height);
            graphicContext.strokeRect(pixelStart, y, pixelSpan, height);
            graphicContext.strokeRect(pixelStart - 1, (y - 1), pixelSpan + 2, height + 2);
            graphicContext.strokeRect(pixelStart - 2, (y - 2), pixelSpan + 4, height + 4);
            if (pixelSpan < 2) {
                graphicContext.strokeRect(pixelStart - 2, (y - 2), pixelSpan + 4, height + 4);
            }
        }
    }

    private ReferenceFrame getReferenceFrame() {
        return frame;
    }
}
