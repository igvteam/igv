/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

import org.broad.igv.feature.SequenceManager;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.util.SOLIDUtils;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class SequenceRenderer implements Renderer {

    static Map<Character, Color> nucleotideColors = new HashMap();

    static {
        //nucleotideColors.put('A', new Color(160, 160, 255));
        //nucleotideColors.put('C', new Color(255, 140, 75));
        //nucleotideColors.put('T', new Color(160, 255, 255));
        //nucleotideColors.put('G', new Color(255, 112, 112));
        //nucleotideColors.put('N', Color.gray);    
        nucleotideColors.put('A', Color.GREEN);
        nucleotideColors.put('a', Color.GREEN);
        nucleotideColors.put('C', Color.BLUE);
        nucleotideColors.put('c', Color.BLUE);
        nucleotideColors.put('T', Color.RED);
        nucleotideColors.put('t', Color.RED);
        nucleotideColors.put('G', new Color(209, 113, 5));
        nucleotideColors.put('g', new Color(209, 113, 5));
        nucleotideColors.put('N', Color.gray);
        nucleotideColors.put('n', Color.gray);
    }

    public void draw(RenderContext context, Rectangle trackRectangle) {
        draw(context, trackRectangle, false);
    }

    // For efficiency
    //@Override
    public void draw(RenderContext context, Rectangle trackRectangle, boolean showColorSpace) {

        double locScale = context.getScale();
        double origin = context.getOrigin();
        String chr = context.getChr();
        String genome = context.getGenomeId();

        int start = Math.max(0, (int) origin - 1);
        int end = (int) (origin + trackRectangle.width * locScale) + 1;
        byte[] seq = SequenceManager.readSequence(genome, chr, start, end);
        byte[] seqCS = null;
        if (showColorSpace) {
            seqCS = SOLIDUtils.convertToColorSpace(seq);
        }

        if (seq != null && seq.length > 0) {
            int hCS = (showColorSpace ? trackRectangle.height / 2 : 0);
            int yBase = hCS + trackRectangle.y + 2;
            int yCS = trackRectangle.y + 2;
            int dY = (showColorSpace ? hCS : trackRectangle.height) - 4;
            int dX = (int) (1.0 / locScale);
            // Create a graphics to use
            Graphics2D g = (Graphics2D) context.getGraphics().create();
            if (dX >= 8) {
                Font f = FontManager.getScalableFont(Font.BOLD, Math.min(dX, 12));
                g.setFont(f);
            }

            // Loop through base pair coordinates
            for (int loc = start; loc < end; loc++) {
                int idx = loc - start;
                int pX0 = (int) ((loc - origin) / locScale);
                char c = (char) seq[idx];
                Color color = nucleotideColors.get(c);
                if (dX >= 8) {
                    if (color == null) {
                        color = Color.black;
                    }
                    g.setColor(color);
                    drawCenteredText(g, new char[]{c}, pX0, yBase + 2, dX, dY - 2);
                    if (showColorSpace) {
                        // draw color space #.  Color space is shifted to be between bases as it represents
                        // two bases.
                        g.setColor(Color.black);
                        String cCS = String.valueOf(seqCS[idx]);
                        drawCenteredText(g, cCS.toCharArray(), pX0 - dX/2, yCS + 2, dX, dY - 2);
                    }
                } else {
                    int bw = dX - 1;
                    if (color != null) {
                        g.setColor(color);
                        g.fillRect(pX0, yBase, bw, dY);
                    }

                }
            }
        }
    }

    private void drawCenteredText(Graphics2D g, char[] chars, int x, int y, int w, int h) {
        // Get measures needed to center the message
        FontMetrics fm = g.getFontMetrics();

        // How many pixels wide is the string
        int msg_width = fm.charsWidth(chars, 0, 1);

        // How far above the baseline can the font go?
        int ascent = fm.getMaxAscent();

        // How far below the baseline?
        int descent = fm.getMaxDescent();

        // Use the string width to find the starting point
        int msgX = x + w / 2 - msg_width / 2;

        // Use the vertical height of this font to find
        // the vertical starting coordinate
        int msgY = y + h / 2 - descent / 2 + ascent / 2;

        g.drawChars(chars, 0, 1, msgX, msgY);

    }

    public void renderBorder(Track track, RenderContext context, Rectangle rect) {


    }

    public void renderAxis(Track track, RenderContext context, Rectangle rect) {


    }

    public void setOverlayMode(boolean overlayMode) {


    }

    public Color getDefaultColor() {
        return Color.BLACK;
    }


}
