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
package org.broad.igv.sam;

import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.SequenceManager;
import org.broad.igv.feature.Strand;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.util.ChromosomeColors;
import org.broad.igv.util.ColorUtilities;

import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class AlignmentRenderer implements FeatureRenderer {

    private static Map<Character, Color> nucleotideColors;
    // Static because all alignment tracks are in color space, or none are
    public static boolean colorSpace;

    private static Color purple = new Color(118, 24, 220);
    private static Color deletionColor = Color.black;
    private static Color skippedColor = new Color(150, 184, 200);
    private static float[] rbgBuffer = new float[3];
    private static float[] colorComps = new float[3];
    private static float[] alignComps = new float[3];
    private static float[] whiteComponents = Color.white.getRGBColorComponents(null);
    private static Color grey2 = new Color(165, 165, 165);
    public static Color grey1 = new Color(200, 200, 200);


    static Font font = FontManager.getScalableFont(10);
    private static Stroke thickStroke = new BasicStroke(2.0f);
    //private static float dash[] = {4.0f, 1.0f};
    //private static Stroke dashedStroke = new BasicStroke(2.0f, BasicStroke.CAP_BUTT,
    //        BasicStroke.JOIN_MITER, 1.0f, dash, 0.0f);

    private static final Color negStrandColor = new Color(110, 145, 225);
    private static final Color posStrandColor = new Color(165, 35, 39);

    private static HashMap<String, Color> readGroupColors = new HashMap();

    PreferenceManager prefs;

    synchronized static private void initColorMap() {
        nucleotideColors = new HashMap();
        nucleotideColors.put('A', Color.GREEN);
        nucleotideColors.put('a', Color.GREEN);
        nucleotideColors.put('C', Color.BLUE);
        nucleotideColors.put('c', Color.BLUE);
        nucleotideColors.put('T', Color.RED);
        nucleotideColors.put('t', Color.RED);
        nucleotideColors.put('G', new Color(209, 113, 5));
        nucleotideColors.put('g', new Color(209, 113, 5));
        nucleotideColors.put('N', Color.gray.brighter());
        nucleotideColors.put('n', Color.gray.brighter());

    }


    /**
     * Constructs ...
     */
    public AlignmentRenderer() {
        this.prefs = PreferenceManager.getInstance();
        if (nucleotideColors == null) {
            initColorMap();
        }

    }

    /**
     * Render a list of alignments in the given rectangle.
     */
    public void renderAlignments(List<Alignment> alignments,
                                 RenderContext context,
                                 Rectangle rect,
                                 AlignmentTrack.RenderOptions renderOptions,
                                 boolean leaveMargin,
                                 Map<String, Color> selectedReadNames) {

        double origin = context.getOrigin();
        double locScale = context.getScale();
        Rectangle screenRect = context.getVisibleRect();

        if ((alignments != null) && (alignments.size() > 0)) {

            //final SAMPreferences prefs = PreferenceManager.getInstance().getSAMPreferences();
            //int insertSizeThreshold = renderOptions.insertSizeThreshold;

            for (Alignment alignment : alignments) {

                // Compute the start and dend of the alignment in pixels
                double pixelStart = ((alignment.getStart() - origin) / locScale);
                double pixelEnd = ((alignment.getEnd() - origin) / locScale);

                // If the any part of the feature fits in the track rectangle draw  it
                if (pixelEnd < rect.x) {
                    continue;
                } else if (pixelStart > rect.getMaxX()) {
                    break;
                }

                Color alignmentColor = getAlignmentColor(alignment, locScale, context.getReferenceFrame().getCenter(), renderOptions);

                Graphics2D g = context.getGraphic2DForColor(alignmentColor);
                g.setFont(font);

                // If the alignment is 3 pixels or less,  draw alignment as posA single block,
                // further detail would not be seen and just add to drawing overhead
                if (pixelEnd - pixelStart < 4) {
                    int w = Math.max(1, (int) (pixelEnd - pixelStart));
                    int h = (int) Math.max(1, rect.getHeight() - 2);
                    int y = (int) (rect.getY() + (rect.getHeight() - h) / 2);
                    g.fillRect((int) pixelStart, y, w, h);
                } else {
                    if (alignment instanceof PairedAlignment) {
                        drawPairedAlignment((PairedAlignment) alignment, rect, g, context, alignmentColor,
                                renderOptions, leaveMargin, selectedReadNames);
                    } else {
                        drawAlignment(alignment, rect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames);
                    }
                }
            }

            // Draw posA border around the center base
            if (locScale < 5 && renderOptions.shadeCenters) {
                // Calculate center lines
                double center = (int) (context.getReferenceFrame().getCenter() - origin);
                int centerLeftP = (int) (center / locScale);
                int centerRightP = (int) ((center + 1) / locScale);
                float transparency = Math.max(0.5f, (float) Math.round(10 * (1 - .75 * locScale)) / 10);
                Graphics2D gBlack = context.getGraphic2DForColor(new Color(0, 0, 0, transparency));
                GraphicUtils.drawDashedLine(gBlack, centerLeftP, rect.y, centerLeftP,
                        rect.y + rect.height);
                if ((centerRightP - centerLeftP > 2)) {
                    GraphicUtils.drawDashedLine(gBlack, centerRightP, rect.y, centerRightP,
                            rect.y + rect.height);
                }
            }
        }
    }

    /**
     * Method for drawing alignments without "blocks" (e.g. DotAlignedAlignment)
     */
    private void drawSimpleAlignment(Alignment alignment,
                                     Rectangle rect,
                                     Graphics2D g,
                                     RenderContext context,
                                     boolean flagUnmappedPair) {
        double origin = context.getOrigin();
        double locScale = context.getScale();
        int x = (int) ((alignment.getStart() - origin) / locScale);
        int length = alignment.getEnd() - alignment.getStart();
        int w = (int) Math.ceil(length / locScale);
        int h = (int) Math.max(1, rect.getHeight() - 2);
        int y = (int) (rect.getY() + (rect.getHeight() - h) / 2);
        int arrowLength = Math.min(5, w / 6);
        int[] xPoly = null;
        int[] yPoly = {y, y, y + h / 2, y + h, y + h};

        // Don't draw off edge of clipping rect
        if (x < rect.x && (x + w) > (rect.x + rect.width)) {
            x = rect.x;
            w = rect.width;
            arrowLength = 0;
        } else if (x < rect.x) {
            int delta = rect.x - x;
            x = rect.x;
            w -= delta;
            if (alignment.isNegativeStrand()) {
                arrowLength = 0;
            }
        } else if ((x + w) > (rect.x + rect.width)) {
            w -= ((x + w) - (rect.x + rect.width));
            if (!alignment.isNegativeStrand()) {
                arrowLength = 0;
            }
        }


        if (alignment.isNegativeStrand()) {
            //     2     1
            //   3
            //     5     5
            xPoly = new int[]{x + w, x, x - arrowLength, x, x + w};
        } else {
            //     1     2
            //             3
            //     5     4
            xPoly = new int[]{x, x + w, x + w + arrowLength, x + w, x};
        }
        g.fillPolygon(xPoly, yPoly, xPoly.length);

        if (flagUnmappedPair && alignment.isPaired() && !alignment.getMate().isMapped()) {
            Graphics2D cRed = context.getGraphic2DForColor(Color.red);
            cRed.drawPolygon(xPoly, yPoly, xPoly.length);
        }
    }

    private void drawPairedAlignment(
            PairedAlignment pair,
            Rectangle rect,
            Graphics2D g,
            RenderContext context,
            Color alignmentColor,
            AlignmentTrack.RenderOptions renderOptions,
            boolean leaveMargin,
            Map<String, Color> selectedReadNames) {

        drawAlignment(pair.firstAlignment, rect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames);
        if (pair.secondAlignment != null) {
            drawAlignment(pair.secondAlignment, rect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames);

            Graphics2D gLine = context.getGraphic2DForColor(grey1);
            double origin = context.getOrigin();
            double locScale = context.getScale();
            int startX = (int) ((pair.firstAlignment.getEnd() - origin) / locScale);
            startX = Math.max(rect.x, startX);

            int endX = (int) ((pair.secondAlignment.getStart() - origin) / locScale);
            endX = Math.min(rect.x + rect.width, endX);

            int h = (int) Math.max(1, rect.getHeight() - (leaveMargin ? 2 : 0));
            int y = (int) (rect.getY()); // + (rect.getHeight() - h) / 2);
            gLine.drawLine(startX, y + h / 2, endX, y + h / 2);

        }

    }

    /**
     * @param alignment
     * @param rect
     * @param g
     * @param context
     * @param alignmentColor
     * @param renderOptions
     * @param leaveMargin
     * @param selectedReadNames
     */
    private void drawAlignment(
            Alignment alignment,
            Rectangle rect,
            Graphics2D g,
            RenderContext context,
            Color alignmentColor,
            AlignmentTrack.RenderOptions renderOptions,
            boolean leaveMargin,
            Map<String, Color> selectedReadNames) {

        double origin = context.getOrigin();
        double locScale = context.getScale();
        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();

        if (blocks == null) {
            drawSimpleAlignment(alignment, rect, g, context, renderOptions.flagUnmappedPairs);
            return;
        }

        AlignmentBlock terminalBlock = alignment.isNegativeStrand() ? blocks[0] : blocks[blocks.length - 1];

        Graphics2D greyGraphics = context.getGraphic2DForColor(new Color(185, 185, 185));

        int lastBlockEnd = Integer.MIN_VALUE;

        int blockNumber = -1;
        char[] gapTypes = alignment.getGapTypes();
        boolean highZoom = locScale < 0.1251;

        for (AlignmentBlock aBlock : alignment.getAlignmentBlocks()) {
            blockNumber++;
            int x = (int) ((aBlock.getStart() - origin) / locScale);
            int w = (int) Math.ceil(aBlock.getBases().length / locScale);
            int h = (int) Math.max(1, rect.getHeight() - (leaveMargin ? 2 : 0));
            int y = (int) (rect.getY()); // + (rect.getHeight() - h) / 2);


            // Create polygon to represent the alignment.
            boolean isZeroQuality = alignment.getMappingQuality() == 0;

            if (highZoom && w > 10) {
                x++;
                w -= 2;
            }

            // If block is out of view skip
            if (x + w >= rect.x && x <= rect.getMaxX()) {
                if (w <= 10 || h <= 10 || aBlock != terminalBlock) {
                    g.fillRect(x, y, w, h);
                    if (isZeroQuality) {
                        greyGraphics.drawRect(x, y, w - 1, h);
                    }

                    if (renderOptions.flagUnmappedPairs && alignment.isPaired() && !alignment.getMate().isMapped()) {
                        Graphics2D cRed = context.getGraphic2DForColor(Color.red);
                        cRed.drawRect(x, y, w, h);
                    }

                    if (selectedReadNames.containsKey(alignment.getReadName())) {
                        Color c = selectedReadNames.get(alignment.getReadName());
                        if (c == null) {
                            c = Color.blue;
                        }
                        Graphics2D cBlue = context.getGraphic2DForColor(c);
                        Stroke s = cBlue.getStroke();
                        cBlue.setStroke(thickStroke);
                        cBlue.drawRect(x, y, w, h);
                        cBlue.setStroke(s);
                    }

                } else {
                    int arrowLength = Math.min(5, w / 6);

                    // Don't draw off edge of clipping rect
                    if (x < rect.x && (x + w) > (rect.x + rect.width)) {
                        x = rect.x;
                        w = rect.width;
                        arrowLength = 0;
                    } else if (x < rect.x) {
                        int delta = rect.x - x;
                        x = rect.x;
                        w -= delta;
                        if (alignment.isNegativeStrand()) {
                            arrowLength = 0;
                        }
                    } else if ((x + w) > (rect.x + rect.width)) {
                        w -= ((x + w) - (rect.x + rect.width));
                        if (!alignment.isNegativeStrand()) {
                            arrowLength = 0;
                        }
                    }

                    int[] xPoly;
                    int[] yPoly = {y, y, y + h / 2, y + h, y + h};

                    if (alignment.isNegativeStrand()) {
                        xPoly = new int[]{x + w, x, x - arrowLength, x, x + w};
                    } else {
                        xPoly = new int[]{x, x + w, x + w + arrowLength, x + w, x};
                    }
                    g.fillPolygon(xPoly, yPoly, xPoly.length);


                    if (isZeroQuality) {
                        greyGraphics.drawPolygon(xPoly, yPoly, xPoly.length);
                    }

                    if (renderOptions.flagUnmappedPairs && alignment.isPaired() && !alignment.getMate().isMapped()) {
                        Graphics2D cRed = context.getGraphic2DForColor(Color.red);
                        cRed.drawPolygon(xPoly, yPoly, xPoly.length);
                    }

                    if (selectedReadNames.containsKey(alignment.getReadName())) {
                        Color c = selectedReadNames.get(alignment.getReadName());
                        if (c == null) {
                            c = Color.blue;
                        }
                        Graphics2D cBlue = context.getGraphic2DForColor(c);
                        Stroke s = cBlue.getStroke();
                        cBlue.setStroke(thickStroke);
                        cBlue.drawPolygon(xPoly, yPoly, xPoly.length);
                        cBlue.setStroke(s);
                    }
                }
            }


            if (locScale < 5) {
                drawBases(context, rect, aBlock, alignmentColor, renderOptions.shadeBases, renderOptions.showAllBases);
            }

            // Draw connecting lines between blocks, if in view
            if (lastBlockEnd > Integer.MIN_VALUE && x > rect.x) {
                Graphics2D gLine;
                Stroke stroke;
                int gapIdx = blockNumber - 1;
                Color gapLineColor = deletionColor;
                if (gapTypes != null && gapIdx < gapTypes.length && gapTypes[gapIdx] == SamAlignment.SKIPPED_REGION) {
                    gLine = context.getGraphic2DForColor(skippedColor);
                    stroke = gLine.getStroke();
                } else {
                    gLine = context.getGraphic2DForColor(gapLineColor);
                    stroke = gLine.getStroke();
                    //gLine.setStroke(dashedStroke);
                    gLine.setStroke(thickStroke);
                }

                int startX = Math.max(rect.x, lastBlockEnd);
                int endX = Math.min(rect.x + rect.width, x);

                gLine.drawLine(startX, y + h / 2, endX, y + h / 2);
                gLine.setStroke(stroke);
            }
            lastBlockEnd = x + w;

            // Next block cannot start before lastBlockEnd.  If its out of view we are done.
            if (lastBlockEnd > rect.getMaxX()) {
                break;
            }

        }

        // Render insertions if locScale ~ 0.25 (base level)
        if (locScale < 0.25) {
            drawInsertions(origin, rect, locScale, alignment, context);
        }
    }

    /**
     * Draw the bases for an alignment block.
     *
     * @param context
     * @param rect
     */
    private void drawBases(RenderContext context,
                           Rectangle rect,
                           AlignmentBlock block,
                           Color alignmentColor,
                           boolean shadeBases,
                           boolean showAllBases) {

        alignmentColor.getRGBColorComponents(alignComps);

        double locScale = context.getScale();
        double origin = context.getOrigin();
        String chr = context.getChr();
        String genome = context.getGenomeId();

        byte[] read = block.getBases();
        boolean isSoftClipped = block.isSoftClipped();

        if ((read != null) && (read.length > 0)) {

            // Compute bounds, get posA graphics to use,  and compute posA font
            int pY = (int) rect.getY();
            int dY = (int) rect.getHeight();
            int dX = (int) Math.max(1, (1.0 / locScale));
            Graphics2D g = (Graphics2D) context.getGraphics().create();
            if (dX >= 8) {
                Font f = FontManager.getScalableFont(Font.BOLD, Math.min(dX, 12));
                g.setFont(f);
            }

            // Get the base qualities, start/end,  and reference sequence

            int start = block.getStart();
            int end = start + read.length;
            byte[] reference = isSoftClipped ? null : SequenceManager.readSequence(genome, chr, start, end);


            // Loop through base pair coordinates
            for (int loc = start; loc < end; loc++) {

                // Index into read array,  just the genomic location offset by
                // the start of this block
                int idx = loc - start;

                // Is this base posA mismatch?  Note '=' means indicates posA match by definition
                // If we do not have posA valid reference we assume posA match.  Soft clipped
                // bases are considered mismatched by definition
                boolean misMatch =
                        isSoftClipped ||
                                (read[idx] != '=' &&
                                        reference != null &&
                                        idx < reference.length &&
                                        reference[idx] != 0 &&
                                        reference[idx] != read[idx]);

                if (misMatch || showAllBases) {
                    char c = (char) read[loc - start];


                    Color color = nucleotideColors.get(c);
                    if (color == null) {
                        color = Color.black;
                    }

                    if (shadeBases) {
                        byte qual = block.getQuality(loc - start);
                        color = getShadedColor(qual, color, prefs);
                    }


                    // If there is room for text draw the character, otherwise
                    // just draw posA rectangle to represent the
                    int pX0 = (int) ((loc - origin) / locScale);

                    // Don't draw out of clipping rect
                    if (pX0 > rect.getMaxX()) {
                        break;
                    } else if (pX0 + dX < rect.getX()) {
                        continue;
                    }

                    if ((dX >= 8) && (dY >= 12)) {
                        g.setColor(color);
                        drawCenteredText(g, new char[]{c}, pX0, pY + 1, dX, dY - 2);
                    } else {

                        int dW = (dX > 4 ? dX - 1 : dX);

                        if (color != null) {
                            g.setColor(color);
                            if (dY < 10) {
                                g.fillRect(pX0, pY, dX, dY);
                            } else {
                                g.fillRect(pX0, pY + 1, dW, dY - 3);
                            }
                        }
                    }
                }

            }
        }
    }

    private Color getShadedColor(byte qual, Color color, PreferenceManager prefs) {
        float alpha = 0;
        int minQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN);
        color.getRGBColorComponents(colorComps);
        if (qual < minQ) {
            alpha = 0.1f;
        } else {
            int maxQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX);
            alpha = Math.max(0.1f, Math.min(1.0f, 0.1f + 0.9f * (qual - minQ) / (maxQ - minQ)));
        }

        // Round alpha to nearest 0.1, for effeciency;
        alpha = ((int) (alpha * 10 + 0.5f)) / 10.0f;
        color = ColorUtilities.getCompositeColor(alignComps, colorComps, alpha);
        return color;
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

    private void drawInsertions(double origin, Rectangle rect, double locScale, Alignment alignment, RenderContext context) {

        Graphics2D gInsertion = context.getGraphic2DForColor(purple);
        AlignmentBlock[] insertions = alignment.getInsertions();
        if (insertions != null) {
            for (AlignmentBlock aBlock : insertions) {
                int x = (int) ((aBlock.getStart() - origin) / locScale);
                int h = (int) Math.max(1, rect.getHeight() - 2);
                int y = (int) (rect.getY() + (rect.getHeight() - h) / 2);

                // Don't draw out of clipping rect
                if (x > rect.getMaxX()) {
                    break;
                } else if (x < rect.getX()) {
                    continue;
                }


                gInsertion.fillRect(x - 2, y, 4, 2);
                gInsertion.fillRect(x - 1, y, 2, h);
                gInsertion.fillRect(x - 2, y + h - 2, 4, 2);
            }
        }
    }

    ContinuousColorScale cs = null;


    private Color getAlignmentColor(Alignment alignment, double locScale,
                                    double center, AlignmentTrack.RenderOptions renderOptions) {

        // Set color used to draw the feature.  Highlight features that intersect the
        // center line.  Also restorePersistentState row "score" if alignment intersects center line

        Color c = grey1;
        switch (renderOptions.colorOption) {
            case INSERT_SIZE:
                boolean isPairedAlignment = alignment instanceof PairedAlignment;
                if (alignment.isPaired() && alignment.getMate().isMapped() || isPairedAlignment) {
                    boolean sameChr = isPairedAlignment ||
                            alignment.getMate().getChr().equals(alignment.getChr());
                    if (sameChr) {
                        int readDistance = Math.abs(alignment.getInferredInsertSize());
                        if (readDistance > 0) {
                            if (readDistance > renderOptions.getMaxInsertSizeThreshold() || readDistance < renderOptions.getMinInsertSizeThreshold()) {
                                c = ChromosomeColors.getColor(alignment.getChr());
                            }
                            //return renderOptions.insertSizeColorScale.getColor(readDistance);
                        }
                    } else {
                        c = ChromosomeColors.getColor(alignment.getMate().getChr());
                        if (c == null) {
                            c = Color.black;
                        }
                    }
                }


                break;
            case PAIR_ORIENTATION:
                c = getOrientationColor(alignment);
                break;
            case READ_STRAND:
                if (alignment.isNegativeStrand()) {
                    c = negStrandColor;
                } else {
                    c = posStrandColor;
                }
                break;
            case FRAGMENT_STRAND:
                if (alignment.getFragmentStrand(1) == Strand.NEGATIVE) {
                    c = negStrandColor;
                } else if (alignment.getFragmentStrand(1) == Strand.POSITIVE) {
                    c = posStrandColor;
                }
                break;
            case READ_GROUP:
                String rg = alignment.getReadGroup();
                if (rg != null) {
                    c = readGroupColors.get(rg);
                    if (c == null) {
                        c = ColorUtilities.randomColor(readGroupColors.size() + 1);
                        readGroupColors.put(rg, c);
                    }
                }
                break;
            case SAMPLE:
                String sample = alignment.getSample();
                if (sample != null) {
                    c = readGroupColors.get(sample);
                    if (c == null) {
                        c = ColorUtilities.randomColor(readGroupColors.size() + 1);
                        readGroupColors.put(sample, c);
                    }
                }
                break;

            default:
                if (renderOptions.shadeCenters && center >= alignment.getStart() && center <= alignment.getEnd()) {
                    if (locScale < 1) {
                        c = grey2;
                    }
                }
        }

        if (alignment.getMappingQuality() == 0) {
            // Maping Q = 0
            float alpha = 0.15f;
            c.getColorComponents(rbgBuffer);
            // Assuming white background TODO -- this should probably be passed in
            return ColorUtilities.getCompositeColor(whiteComponents, rbgBuffer, alpha);
        }

        return c;
    }


    /**
     * Illumin scheme -- todo, something for Solid
     *
     * @return
     */

    private Color getOrientationColor(Alignment alignment) {

        Color c = null;
        if (!alignment.isProperPair() && !alignment.isSmallInsert()) {
            if (alignment.getAttribute("CS") != null) {
                c = getSolidScheme().get(alignment.getPairOrientation());
            } else {
                c = getIlluminaScheme().get(alignment.getPairOrientation());
            }
        }

        return c == null ? grey1 : c;

    }

    private static final Color LR_COLOR = grey1; // "Normal" alignment color
    private static final Color RL_COLOR = new Color(0, 150, 0);
    private static final Color RR_COLOR = new Color(0, 0, 150);
    private static final Color LL_COLOR = new Color(0, 150, 150);

    Map<String, Color> illuminaScheme;
    Map<String, Color> solidScheme;

    private Map<String, Color> getIlluminaScheme() {
        if (illuminaScheme == null) {
            illuminaScheme = new HashMap();
            //LR
            illuminaScheme.put("F1R2", LR_COLOR);
            illuminaScheme.put("F2R1", LR_COLOR);
            //LL
            illuminaScheme.put("F1F2", LL_COLOR);
            illuminaScheme.put("F2F1", LL_COLOR);
            //RR
            illuminaScheme.put("R1R2", RR_COLOR);
            illuminaScheme.put("R2R1", RR_COLOR);
            //RL
            illuminaScheme.put("R1F2", RL_COLOR);
            illuminaScheme.put("R2F1", RL_COLOR);
        }
        return illuminaScheme;
    }


    private Map<String, Color> getSolidScheme() {
        if (solidScheme == null) {
            solidScheme = new HashMap();
            //LR
            solidScheme.put("F1F2", LR_COLOR);
            solidScheme.put("R2R1", LR_COLOR);
            //LL -- switched with RR color per Bob's instructions
            solidScheme.put("F1R2", RR_COLOR);
            solidScheme.put("R2F1", RR_COLOR);
            //RR
            solidScheme.put("R1F2", LL_COLOR);
            solidScheme.put("F2R1", LL_COLOR);
            //RL
            solidScheme.put("R1R2", RL_COLOR);
            solidScheme.put("F2F1", RL_COLOR);
        }
        return solidScheme;
    }


    public static Map<Character, Color> getNucleotideColors() {
        if (nucleotideColors == null) {
            initColorMap();
        }
        return nucleotideColors;
    }
}
