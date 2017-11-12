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
 * CytobandRenderer.java
 *
 * Created on September 19, 2007, 5:36 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.ui.javafx.renderer;

import com.sun.javafx.tk.FontMetrics;
import javafx.geometry.Rectangle2D;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.paint.Color;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CytobandRenderer {

    boolean drawLabels = true;
    static final public double CYTOBAND_Y_OFFSET = 5.0;
    static final public double LABEL_OFFSET = 25.0;
    private static Map<Integer, Color> stainColors = new HashMap<Integer, Color>();
    private FontMetrics fontMetrics;

    public CytobandRenderer(FontMetrics fontMetrics) {
        this.fontMetrics = fontMetrics;
    }

    public void draw(List<Cytoband> data, GraphicsContext graphicsContext, Rectangle2D graphicRect, ReferenceFrame frame) {

        if (data.size() > 0) {

            // TODO: not dealing with image export yet

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
            for (Cytoband cytoband : cytobands) {
                double s = sc * cytoband.getStart();
                double e = sc * cytoband.getEnd();
                double stringWidth = fontMetrics.computeStringWidth(cytoband.getName());
                double x = s + (e - s - stringWidth) / 2;
                if (x > (prevEnd + minSpacing)) {
                    graphicsContext.strokeText(cytoband.getName(), x, adjustedY);
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
