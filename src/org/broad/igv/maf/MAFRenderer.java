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
package org.broad.igv.maf;

import org.broad.igv.maf.MAFTile.Gap;
import org.broad.igv.maf.MAFTile.MASequence;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;

import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class MAFRenderer {

    static Map<Character, Color> nucleotideColors = new HashMap();


    static {

        nucleotideColors.put('A', Color.GREEN);
        nucleotideColors.put('a', Color.GREEN);
        nucleotideColors.put('C', Color.BLUE);
        nucleotideColors.put('c', Color.BLUE);
        nucleotideColors.put('T', Color.RED);
        nucleotideColors.put('t', Color.RED);
        nucleotideColors.put('G', new Color(242, 182, 65));
        nucleotideColors.put('g', new Color(242, 182, 65));
        nucleotideColors.put('N', Color.gray);
        nucleotideColors.put('n', Color.gray);
    }

    /**
     * Constructs ...
     */
    public MAFRenderer() {
    }

    public void renderGaps(List<Gap> gaps, RenderContext context, Rectangle rect) {
        double origin = context.getOrigin();
        double locScale = context.getScale();

        if (locScale > 1) {
            return;
        }

        for (Gap gap : gaps) {

            int pixelPosition = (int) ((gap.getPosition() - origin) / locScale);

            Graphics2D g = context.getGraphic2DForColor(Color.BLACK);

            g.drawLine(pixelPosition, rect.y + rect.height - 5, pixelPosition,
                    rect.y + rect.height);

            Rectangle textRect = new Rectangle(rect);
            textRect.x = pixelPosition - 6;
            textRect.width = 12;
            textRect.height = rect.height - 7;
            GraphicUtils.drawCenteredText(String.valueOf(gap.getCount()),
                    textRect, g);
        }
    }

    /**
     * Method description
     *
     * @param context
     * @param trackRectangle
     * @param track
     */
    public void renderAligment(
            MASequence alignedSequence,
            MASequence reference,
            List<Gap> gaps,
            RenderContext context,
            Rectangle trackRectangle,
            Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        if (locScale > 1) {
            return;
        }

// Note -- don't cast these to an int until the range is checked.
// could get an overflow.
        double pixelStart = ((alignedSequence.getStart() - origin) / locScale);
        double pixelEnd = ((alignedSequence.getEnd() - origin) / locScale);

        // If the any part of the alignedSequence fits in the track rectangle draw it
        if ((pixelEnd >= trackRectangle.getX()) && (pixelStart <= trackRectangle.getMaxX())) {


            int w = Math.max(1, (int) (pixelEnd - pixelStart));
            int h = (int) Math.max(1, trackRectangle.getHeight() - 2);
            int y = (int) (trackRectangle.getY() + (trackRectangle.getHeight() - h) / 2);
            Rectangle block = new Rectangle((int) pixelStart, y, w, h);

            if (locScale < 1) {
                drawBases(context, block, alignedSequence, reference, gaps);
            }

        }

    }

    /**
     * Method description
     *
     * @param context
     * @param rect
     * @param alignment
     */
    private void drawBases(RenderContext context, Rectangle rect,
                           MASequence alignment, MASequence ref, List<Gap> gaps) {

        double locScale = context.getScale();
        double origin = context.getOrigin();

        {

            int pY = (int) rect.getY();
            int dY = (int) rect.getHeight();
            int dX = (int) (1.0 / locScale);

            // Create a graphics to use
            Graphics2D g = (Graphics2D) context.getGraphics().create();
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            if (dX >= 8) {
                Font f = FontManager.getFont(Font.BOLD, Math.min(dX, 12));
                g.setFont(f);
            }

            // Loop through base pair coordinates
            int windowStart = (int) origin - 1;
            int windowEnd = (int) context.getEndLocation() + 1;
            int start = Math.max(windowStart, alignment.getStart());
            int end = Math.min(windowEnd, alignment.getEnd());

            for (int loc = start; loc < end; loc++) {

                int pX0 = (int) ((loc - origin) / locScale);

                char c = alignment.getGapAdjustedBase(loc);
                char refBase = ref.getGapAdjustedBase(loc);
                if (c == 0 || refBase == 0) {
                    continue;
                }

                boolean misMatch = Character.toUpperCase(c) != Character.toUpperCase(refBase);

                char charToDraw = misMatch || ref == alignment ? c : '.';

                Color color = nucleotideColors.get(charToDraw);

                if ((dX >= 8) && (dY >= 12) || charToDraw == '.') {

                    // Graphics2D gBackground = context.getGraphic2DForColor(background);
                    // gBackground.fillRect(pX0, pY, dX, dY);
                    if (charToDraw == '.') {
                        color = Color.LIGHT_GRAY;
                    } else if (color == null) {
                        color = Color.black;
                    }

                    g.setColor(color);
                    drawCenteredText(g, new char[]{charToDraw}, pX0, pY + 2, dX, dY - 2);
                } else {
                    if (color != null) {
                        g.setColor(color);
                        g.fillRect(pX0, pY, dX - 1, dY);
                    }

                }
            }

            // Check for insertion
            if (gaps != null) {
                Graphics2D gapG = context.getGraphic2DForColor(Color.black);
                for (Gap gap : gaps) {
                    for (int idx = gap.startIdx; idx < gap.endIdx; idx++) {
                        if (alignment.bases.charAt(idx) != '-') {
                            int pX0 = (int) ((gap.getPosition() - origin) / locScale);
                            gapG.drawLine(pX0, pY, pX0, pY + dY);
                            break;
                        }
                    }
                }

            }
        }
    }

    private void drawCenteredText(Graphics2D g, char[] chars, int x, int y,
                                  int w, int h) {

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


}

