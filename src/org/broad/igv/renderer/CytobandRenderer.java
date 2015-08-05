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
package org.broad.igv.renderer;

import org.broad.igv.feature.Cytoband;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class CytobandRenderer {

    boolean drawLabels = true;
    static final public int CYTOBAND_Y_OFFSET = 5;
    static final public int LABEL_OFFSET = 25;
    private static Map<Integer, Color> stainColors = new HashMap<Integer, Color>();

    public void draw(List<Cytoband> data, Graphics g2D, Rectangle graphicRect, ReferenceFrame frame) {

        if (data.size() > 0) {

            // If we are in the process of exporting to an image file
            // we need to write out the cytoband locus on the image
            if (IGV.getInstance().isExportingSnapshot() || FrameManager.isGeneListMode()) {
                String locus = frame.getChrName();
                if (locus != null) {
                    Graphics g2 = g2D.create();
                    g2.setFont(FontManager.getFont(Font.BOLD, 11));
                    g2.drawString(locus, 3, 11);
                    g2.dispose();
                }

            }

            // Draw Cytoband
            drawBands(data, g2D, graphicRect, frame.getMaxCoordinate());

            // Draw Cytoband Labels
            if (drawLabels && !FrameManager.isGeneListMode()) {
                drawLabels(g2D, graphicRect, data, frame.getMaxCoordinate());
            }
        }
    }

    public void drawBands(List<Cytoband> data, Graphics g2D, Rectangle graphicRect, int chromoLength) {

        int[] xC = new int[3];
        int[] yC = new int[3];

        double scale = graphicRect.getWidth() / chromoLength;

        int lastPX = -1;
        for (Cytoband cytoband : data) {

            int start = (int) (graphicRect.getX() + scale * cytoband.getStart());
            int end = (int) (graphicRect.getX() + scale * cytoband.getEnd());
            if (end > lastPX) {

                int y = (int) graphicRect.getY() + CYTOBAND_Y_OFFSET;
                int height = (int) graphicRect.getHeight();

                if (cytoband.getType() == 'c') { // centermere: "acen"

                    int center = (y + height / 2);
                    if (cytoband.getName().startsWith("p")) {
                        xC[0] = start;
                        yC[0] = (int) graphicRect.getMaxY() + CYTOBAND_Y_OFFSET;
                        xC[1] = start;
                        yC[1] = y;
                        xC[2] = end;
                        yC[2] = center;
                    } else {
                        xC[0] = end;
                        yC[0] = (int) graphicRect.getMaxY() + CYTOBAND_Y_OFFSET;
                        xC[1] = end;
                        yC[1] = y;
                        xC[2] = start;
                        yC[2] = center;
                    }
                    g2D.setColor(Color.RED.darker());
                    g2D.fillPolygon(xC, yC, 3);
                } else {

                    g2D.setColor(getCytobandColor(cytoband));
                    g2D.fillRect(start, y, (end - start), height);
                    g2D.setColor(Color.BLACK);
                    g2D.drawRect(start, y, (end - start), height);
                }
            }
            lastPX = end;
        }
    }

    private void drawLabels(final Graphics g, Rectangle graphicRect, List<Cytoband> cytobands, int chromoLength) {

        double width = graphicRect.getWidth();
        int y = (int) graphicRect.getY() + LABEL_OFFSET;

        // Draw labels
        g.setColor(Color.BLACK);
        FontMetrics fm = g.getFontMetrics();
        int minSpacing = 10;
        int prevEnd = 0;
        double sc = width / chromoLength;
        int adjustedY = y;
        if (cytobands != null) {
            for (Cytoband cytoband : cytobands) {
                int s = (int) (sc * cytoband.getStart());
                int e = (int) (sc * cytoband.getEnd());
                int stringWidth = (int) fm.getStringBounds(cytoband.getName(), g).getWidth();
                int x = (int) (s + (e - s - stringWidth) / 2);
                if (x > (prevEnd + minSpacing)) {
                    g.drawString(cytoband.getName(), x, adjustedY);
                    prevEnd = x + stringWidth;
                }
            }
        }
    }

    public void applyTextTranslationAndRotation(Graphics2D graphics2, int x, int y) {

        AffineTransform transform = new AffineTransform();
        transform.translate(x, y);
        transform.rotate(-Math.PI / 2);
        graphics2.transform(transform);
    }

    private static Color getCytobandColor(Cytoband data) {
        if (data.getType() == 'p') { // positive: "gpos100"
            int stain = data.getStain(); // + 4;

            int shade = (int) (255 - stain / 100.0 * 255);
            Color c = stainColors.get(shade);
            if (c == null) {
                c = new Color(shade, shade, shade);
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
