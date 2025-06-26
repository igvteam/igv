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

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.SnapshotUtilities;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class CytobandRenderer {

    private final boolean darkMode;
    boolean drawLabels = true;
    static final public int CYTOBAND_Y_OFFSET = 5;
    static final public int LABEL_OFFSET = 25;
    private int fontHeight = 10;
    private int bandHeight = 10;
    private String fontFamilyName = "Lucida Sans";

    private static Map<Integer, Color> stainColors = new HashMap<Integer, Color>();

    public CytobandRenderer(boolean darkMode) {
        this.darkMode = darkMode;

//        if (darkMode) {
//            fontFamilyName = "Lucida Sans Typewriter";
//        }
//        // Initialize stain colors
//        for (int i = 0; i <= 100; i += 5) {
//            int shade = (int) (255 - i / 100.0 * 255);
//            stainColors.put(i, new Color(shade, shade, shade));
//        }
    }

    public void drawIdeogram(List<Cytoband> data, Graphics g2D, Rectangle graphicRect, ReferenceFrame frame) {

        if (data.size() > 0) {
            Graphics g = g2D.create();
            try {

                FontManager.getFont(fontHeight);
                g.setFont(new Font(fontFamilyName, Font.BOLD, fontHeight));

                // If we are in the process of exporting to an image file
                // we need to write out the cytoband locus on the image
                if (SnapshotUtilities.snapshotInProgress) {
                    String locus = frame.getFormattedLocusString();
                    if (locus != null) {
                        Graphics g2 = g2D.create();
                        if(darkMode) {
                            g2.setColor(Color.WHITE);
                        }
                        g2.setFont(FontManager.getFont(Font.BOLD, 11));
                        g2.drawString(locus, 3, 11);
                        g2.dispose();
                    }
                }

                // Draw Cytoband
                drawBand(data, g, graphicRect, 0, frame.getMaxCoordinate());

                // Draw Cytoband Labels
                if (drawLabels && !FrameManager.isGeneListMode()) {
                    drawLabels(g, graphicRect, data, 0, frame.getMaxCoordinate());
                }
            } finally {
                g.dispose();
            }
        }
    }

    public void drawTrack(List<Cytoband> data, Graphics g2D, Rectangle graphicRect, ReferenceFrame frame) {

        if (data.size() > 0) {
            Graphics g = g2D.create();
            try {

                FontManager.getFont(fontHeight);
                g.setFont(new Font(fontFamilyName, Font.BOLD, fontHeight));
                Rectangle cytoRect = new Rectangle(0, 0, graphicRect.width, 10);

                double start = frame.getOrigin();
                double end = frame.getEnd();
                drawBand(data, g, cytoRect, start, end);
                drawLabels(g, graphicRect, data, start, end);
            } finally {
                g.dispose();
            }
        }
    }

    private void drawBand(List<Cytoband> data, Graphics g2D, Rectangle graphicRect, double start, double end) {

        int[] xC = new int[3];
        int[] yC = new int[3];

        double scale = graphicRect.getWidth() / (end - start);

        int lastPX = -1;
        for (Cytoband cytoband : data) {
            if(cytoband.getEnd() < start) continue;
            if(cytoband.getStart() > end) break;
            int s = (int) (graphicRect.getX() + scale * (cytoband.getStart() - start));
            int e = (int) (graphicRect.getX() + scale * (cytoband.getEnd() - start));
            if (e > lastPX) {

                int y = (int) graphicRect.getY() + CYTOBAND_Y_OFFSET;
                int height = (int) graphicRect.getHeight();

                if (cytoband.getType() == 'c') { // centermere: "acen"

                    int center = (y + height / 2);
                    if (cytoband.getName().startsWith("p")) {
                        xC[0] = s;
                        yC[0] = (int) graphicRect.getMaxY() + CYTOBAND_Y_OFFSET;
                        xC[1] = s;
                        yC[1] = y;
                        xC[2] = e;
                        yC[2] = center;
                    } else {
                        xC[0] = e;
                        yC[0] = (int) graphicRect.getMaxY() + CYTOBAND_Y_OFFSET;
                        xC[1] = e;
                        yC[1] = y;
                        xC[2] = s;
                        yC[2] = center;
                    }
                    g2D.setColor(Color.RED.darker());
                    g2D.fillPolygon(xC, yC, 3);
                } else {

                    g2D.setColor(getCytobandColor(cytoband));
                    g2D.fillRect(s, y, (e - s), height);
                    g2D.setColor(Color.BLACK);
                    g2D.drawRect(s, y, (e - s), height);
                }
            }
            lastPX = e;
        }
    }

    private void drawLabels(final Graphics g, Rectangle graphicRect, List<Cytoband> cytobands, double start,  double end) {

        double width = graphicRect.getWidth();
        int y = (int) graphicRect.getY() + LABEL_OFFSET;

        // Draw labels
        g.setColor(Color.BLACK);
        FontMetrics fm = g.getFontMetrics();
        int minSpacing = 10;
        int prevEnd = 0;
        double sc = width / (end - start);
        int adjustedY = y;
        if (cytobands != null) {
            for (Cytoband cytoband : cytobands) {
                if(cytoband.getEnd() < start) continue;
                if(cytoband.getStart() > end) break;
                int s = (int) (sc * (cytoband.getStart() - start));
                int e = (int) (sc * (cytoband.getEnd() - start));
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
