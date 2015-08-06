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

package org.broad.igv.maf;

import org.broad.igv.PreferenceManager;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.sam.AlignmentRenderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;

import java.awt.*;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class MAFRenderer {

//    static Map<Character, Color> nucleotideColors = new HashMap();
//
//
//    static {
//
//        nucleotideColors.put('A', Color.GREEN);
//        nucleotideColors.put('a', Color.GREEN);
//        nucleotideColors.put('C', Color.BLUE);
//        nucleotideColors.put('c', Color.BLUE);
//        nucleotideColors.put('T', Color.RED);
//        nucleotideColors.put('t', Color.RED);
//        nucleotideColors.put('G', new Color(242, 182, 65));
//        nucleotideColors.put('g', new Color(242, 182, 65));
//        nucleotideColors.put('N', Color.gray);
//        nucleotideColors.put('n', Color.gray);
//    }

    /**
     * Constructs ...
     */
    public MAFRenderer() {
    }

    public void renderGaps(List<MultipleAlignmentBlock.Gap> gaps, RenderContext context, Rectangle rect) {
        double origin = context.getOrigin();
        double locScale = context.getScale();

        if (locScale > 1) {
            return;
        }

        for (MultipleAlignmentBlock.Gap gap : gaps) {

            int pixelPosition = (int) ((gap.position - origin) / locScale);

            Graphics2D g = context.getGraphic2DForColor(Color.BLACK);

            g.drawLine(pixelPosition, rect.y + rect.height - 5, pixelPosition,
                    rect.y + rect.height);

            Rectangle textRect = new Rectangle(rect);
            textRect.x = pixelPosition - 6;
            textRect.width = 12;
            textRect.height = rect.height - 7;
            GraphicUtils.drawCenteredText(String.valueOf(gap.size),
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
    public void renderSequence(
            MultipleAlignmentBlock multipleAlignment,
            MultipleAlignmentBlock.Sequence alignedSequence,
            MultipleAlignmentBlock.Sequence reference,
            List<MultipleAlignmentBlock.Gap> gaps,
            RenderContext context,
            Rectangle trackRectangle,
            Track track) {

        Map<Character, Color> nucleotideColors = SequenceRenderer.getNucleotideColors();

        double origin = context.getOrigin();
        double locScale = context.getScale();

        if (locScale > 1) {
            return;
        }

        double pixelStart = ((multipleAlignment.getStart() - origin) / locScale);


        int w = Math.max(1, (int) (trackRectangle.width));
        int h = (int) Math.max(1, trackRectangle.getHeight() - 2);
        int y = (int) (trackRectangle.getY() + (trackRectangle.getHeight() - h) / 2);
        Rectangle rect = new Rectangle((int) pixelStart, y, w, h);

        if (locScale < 1) {

            int pY = (int) rect.getY();
            int dY = (int) rect.getHeight();
            int dX = (int) (1.0 / locScale);

            // Create a graphics to use
            Graphics2D g = (Graphics2D) context.getGraphics().create();
            if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
                g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            }

            if (dX >= 8) {
                Font f = FontManager.getFont(Font.BOLD, Math.min(dX, 12));
                g.setFont(f);
            }

            // Loop through base pair coordinates
            int windowStart = (int) origin - 1;
            int windowEnd = (int) context.getEndLocation() + 1;
            int start = Math.max(windowStart, multipleAlignment.getStart());
            int end = Math.min(windowEnd, multipleAlignment.getEnd());

            byte[] alignmentBytes = alignedSequence.getText().getBytes();
            byte[] refBytes = reference.getText().getBytes();

            for (int loc = start; loc < end; loc++) {

                int pX0 = (int) ((loc - origin) / locScale);

                int idx = multipleAlignment.getGapAdjustedIndex(loc);

                char c = (char) alignmentBytes[idx];
                char refBase = (char) refBytes[idx];

                boolean misMatch = Character.toUpperCase(c) != Character.toUpperCase(refBase);

                char charToDraw = misMatch || reference == alignedSequence ? c : '.';

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
                for (MultipleAlignmentBlock.Gap gap : gaps) {
                    for (int idx = gap.startIdx; idx < gap.startIdx + gap.size; idx++) {
                        if (alignmentBytes[idx] != '-') {
                            int pX0 = (int) ((gap.position - origin) / locScale);
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

