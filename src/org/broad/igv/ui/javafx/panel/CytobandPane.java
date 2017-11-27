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

import javafx.geometry.Rectangle2D;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.ui.javafx.FontMetrics;
import org.broad.igv.ui.javafx.ResizableCanvas;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

//~--- non-JDK imports --------------------------------------------------------

public class CytobandPane extends ResizableCanvas {

    private static Logger log = Logger.getLogger(CytobandPane.class);
    private static double fontHeight = 10.0;
    private static final double bandHeight = 10.0;
    private static String fontFamilyName = "Lucida Sans";
    static final public double CYTOBAND_Y_OFFSET = 5.0;
    static final public double LABEL_OFFSET = 25.0;

    private boolean isDragging = false;

    private double viewOrigin;
    private double viewEnd;

    private double cytobandScale;
    private ReferenceFrame frame;
    private List<Cytoband> currentCytobands;

    private boolean drawLabels = true;
    private static Map<Integer, Color> stainColors = new HashMap<Integer, Color>();

    public CytobandPane(ReferenceFrame frame) {
        this(frame, true);
    }

    public CytobandPane(ReferenceFrame frame, boolean mouseable) {

        this.frame = frame;
        viewOrigin = frame.getOrigin();
        viewEnd = frame.getEnd();

        Font font = Font.font(fontFamilyName, FontWeight.BOLD, fontHeight);
        getCanvas().getGraphicsContext2D().setFont(font);
        setMinHeight(40);
        setMaxHeight(40);
        setPrefHeight(40);

        // Re-render on change of width/height or chromosome.  Note that the height is fixed and
        // so that listener should never execute.  However, leaving this in place as a pattern
        // and also because it may not be fixed in the long run.
        frame.chromosomeNameProperty().addListener((observable, oldValue, newValue) -> render());
        this.prefWidthProperty().addListener((observable, oldValue, newValue) -> render());
        this.prefHeightProperty().addListener((observable, oldValue, newValue) -> render());
        
        render();
    }

    private void render() {
        log.info("rendering " + frame.getChrName());
        log.info("cytoband HW: " + getWidth() + ":" + getHeight());
        log.info("cytoband pHW: " + getPrefWidth() + ":" + getPrefHeight());

        Canvas canvas = getCanvas();
        GraphicsContext graphicContext = canvas.getGraphicsContext2D();
        graphicContext.clearRect(0.0, 0.0, canvas.getWidth(), canvas.getHeight());

        log.info("cytoband canvas HW: " + canvas.getWidth() + ":" + canvas.getHeight());
        
        if (frame.getChrName().equals(Globals.CHR_ALL)) {
            log.info("Frame set at ALL; exiting layoutChildren");
            return;
        }

        Chromosome chromosome = getReferenceFrame().getChromosome();
        if (chromosome == null) {
            log.info("No chromosome; exiting layoutChildren");
            return;
        }

        currentCytobands = chromosome.getCytobands();
        if (currentCytobands == null) {
            log.info("No cytobands; exiting layoutChildren");
            return;
        }

        // The original Swing code gets the width from the Frame, but it seems more straightforward 
        // to get it from the Pane itself.  Do it this way for now, but possibly change later.
        //int dataPanelWidth = frame.getWidthInPixels();
        double dataPanelWidth = this.getPrefWidth();
        Rectangle2D cytoRect = new Rectangle2D(0.0, 10.0, dataPanelWidth, bandHeight);

        log.info("drawing cytoband");
        draw(currentCytobands, graphicContext, cytoRect, frame);
        log.info("done drawing cytoband");

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

            double y = cytoRect.getMinY() + CYTOBAND_Y_OFFSET;
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


    public void draw(List<Cytoband> data, GraphicsContext graphicsContext, Rectangle2D graphicRect, ReferenceFrame frame) {

        if (data.size() > 0) {

            // TODO: image export.  Not dealing with this yet

            // Draw Cytoband
            drawBands(data, graphicsContext, graphicRect, frame.getMaxCoordinate());

            // Draw Cytoband Labels
            if (drawLabels && !FrameManager.isGeneListMode()) {
                drawLabels(graphicsContext, graphicRect, data, frame.getMaxCoordinate());
            }
        }
    }

    private void drawBands(List<Cytoband> data, GraphicsContext graphicsContext, Rectangle2D graphicRect, int chromoLength) {

        double[] xC = new double[3];
        double[] yC = new double[3];

        double scale = graphicRect.getWidth() / chromoLength;

        double lastPX = -1;
        for (Cytoband cytoband : data) {

            double start = graphicRect.getMinX() + scale * cytoband.getStart();
            double end = graphicRect.getMinX() + scale * cytoband.getEnd();
            if (end > lastPX) {

                double y = graphicRect.getMinY() + CYTOBAND_Y_OFFSET;
                double height = graphicRect.getHeight();

                if (cytoband.getType() == 'c') { // centermere: "acen"

                    double center = y + height / 2;
                    if (cytoband.getName().startsWith("p")) {
                        xC[0] = start;
                        yC[0] = graphicRect.getMaxY() + CYTOBAND_Y_OFFSET;
                        xC[1] = start;
                        yC[1] = y;
                        xC[2] = end;
                        yC[2] = center;
                    } else {
                        xC[0] = end;
                        yC[0] = graphicRect.getMaxY() + CYTOBAND_Y_OFFSET;
                        xC[1] = end;
                        yC[1] = y;
                        xC[2] = start;
                        yC[2] = center;
                    }
                    graphicsContext.setFill(Color.RED.darker());
                    graphicsContext.fillPolygon(xC, yC, 3);
                } else {

                    graphicsContext.setFill(getCytobandColor(cytoband));
                    graphicsContext.fillRect(start, y, (end - start), height);
                    graphicsContext.setStroke(Color.BLACK);
                    graphicsContext.strokeRect(start, y, (end - start), height);
                }
            }
            lastPX = end;
        }
    }

    private void drawLabels(final GraphicsContext graphicsContext, Rectangle2D graphicRect, List<Cytoband> cytobands, int chromoLength) {

        double width = graphicRect.getWidth();
        double y = graphicRect.getMinY() + LABEL_OFFSET;

        // Draw labels
        graphicsContext.setStroke(Color.BLACK);
        double minSpacing = 10.0;
        double prevEnd = 0.0;
        double sc = width / chromoLength;
        double adjustedY = y;
        if (cytobands != null) {
            Text sizer = new Text("");
            sizer.setFont(graphicsContext.getFont());
            for (Cytoband cytoband : cytobands) {
                double s = sc * cytoband.getStart();
                double e = sc * cytoband.getEnd();
                double stringWidth = FontMetrics.getTextWidthInFont(cytoband.getName(), sizer);
                double x = s + (e - s - stringWidth) / 2;
                if (x > (prevEnd + minSpacing)) {
                    graphicsContext.setFill(Color.BLACK);
                    graphicsContext.fillText(cytoband.getName(), x, adjustedY);
                    prevEnd = x + stringWidth;
                }
            }
        }
    }

    private static Color getCytobandColor(Cytoband data) {
        if (data.getType() == 'p') {
            int stain = data.getStain();

            int shade = (int) (255 - stain / 100.0 * 255);
            Color c = stainColors.get(shade);
            if (c == null) {
                c = Color.rgb(shade, shade, shade);
                stainColors.put(shade, c);
            }

            return c;

        } else if (data.getType() == 'c') { // centermere: "acen"
            return Color.PINK;

        } else {
            return Color.WHITE;
        }
    }
}
