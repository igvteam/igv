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
package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.sam.AlignmentTrack.ColorOption;
import org.broad.igv.sam.AlignmentTrack.RenderOptions;
import org.broad.igv.sam.BisulfiteBaseInfo.DisplayStatus;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ChromosomeColors;

import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class AlignmentRenderer implements FeatureRenderer {

    private static Logger log = Logger.getLogger(AlignmentRenderer.class);

    private static Map<Character, Color> nucleotideColors;

    private static Color smallISizeColor = new Color(0, 0, 150);
    private static Color largeISizeColor = new Color(150, 0, 0);
    private static Color purple = new Color(118, 24, 220);
    private static Color deletionColor = Color.black;
    private static Color skippedColor = new Color(150, 184, 200);
    public static Color grey1 = new Color(200, 200, 200);

    private static Stroke thickStroke = new BasicStroke(2.0f);

    public static final Color negStrandColor = new Color(150, 150, 200);
    public static final Color posStrandColor = new Color(200, 150, 150);

    private static HashMap<String, Color> readGroupColors = new HashMap();

    private static final Color LR_COLOR = grey1; // "Normal" alignment color
    private static final Color RL_COLOR = new Color(0, 150, 0);
    private static final Color RR_COLOR = new Color(0, 0, 150);
    private static final Color LL_COLOR = new Color(0, 150, 150);
    private static final Color OUTLINE_COLOR = new Color(185, 185, 185);
    public static final Color GROUP_DIVIDER_COLOR = new Color(200, 200, 200);

    static Map<String, Color> frOrientationColors;
    static Map<String, Color> ffOrientationColors;
    static Map<String, Color> rfOrientationColors;

    // Bisulfite constants
    private static final Color bisulfiteColorFw1 = new Color(195, 195, 195); // A bit lighter than normal LR_COLOR
    private static final Color bisulfiteColorRev1 = new Color(195, 222, 195); // A bit lighter than normal LR_COLOR
    private static final Color bisulfiteColorRev2 = new Color(210, 210, 210); // A bit lighter than normal LR_COLOR
    private static final Color bisulfiteColorFw2 = new Color(210, 222, 210); // A bit lighter than normal LR_COLOR
    private static final Color nomeseqColor = new Color(195, 195, 195); // A bit lighter than normal LR_COLOR

    PreferenceManager prefs;

    static {
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


        // fr Orienations (e.g. Illumina paired-end libraries)
        frOrientationColors = new HashMap();
        //LR
        frOrientationColors.put("F1R2", LR_COLOR);
        frOrientationColors.put("F2R1", LR_COLOR);
        frOrientationColors.put("F R ", LR_COLOR);
        frOrientationColors.put("FR", LR_COLOR);
        //LL
        frOrientationColors.put("F1F2", LL_COLOR);
        frOrientationColors.put("F2F1", LL_COLOR);
        frOrientationColors.put("F F ", LL_COLOR);
        frOrientationColors.put("FF", LL_COLOR);
        //RR
        frOrientationColors.put("R1R2", RR_COLOR);
        frOrientationColors.put("R2R1", RR_COLOR);
        frOrientationColors.put("R R ", RR_COLOR);
        frOrientationColors.put("RR", RR_COLOR);
        //RL
        frOrientationColors.put("R1F2", RL_COLOR);
        frOrientationColors.put("R2F1", RL_COLOR);
        frOrientationColors.put("R F ", RL_COLOR);
        frOrientationColors.put("RF", RL_COLOR);

        // rf orienation  (e.g. Illumina mate-pair libraries)
        rfOrientationColors = new HashMap();
        //LR
        rfOrientationColors.put("R1F2", LR_COLOR);
        rfOrientationColors.put("R2F1", LR_COLOR);
        rfOrientationColors.put("R F ", LR_COLOR);
        rfOrientationColors.put("RF", LR_COLOR);
        //LL
        rfOrientationColors.put("R1R2", LL_COLOR);
        rfOrientationColors.put("R2R1", LL_COLOR);
        rfOrientationColors.put("R R ", LL_COLOR);
        rfOrientationColors.put("RR ", LL_COLOR);

        rfOrientationColors.put("F1F2", RR_COLOR);
        rfOrientationColors.put("F2F1", RR_COLOR);
        rfOrientationColors.put("F F ", RR_COLOR);
        rfOrientationColors.put("FF", RR_COLOR);
        //RL
        rfOrientationColors.put("F1R2", RL_COLOR);
        rfOrientationColors.put("F2R1", RL_COLOR);
        rfOrientationColors.put("F R ", RL_COLOR);
        rfOrientationColors.put("FR", RL_COLOR);


        // ff orienation  (e.g. SOLID libraries)
        ffOrientationColors = new HashMap();
        //LR
        ffOrientationColors.put("F1F2", LR_COLOR);
        ffOrientationColors.put("R2R1", LR_COLOR);
        //LL -- switched with RR color per Bob's instructions
        ffOrientationColors.put("F1R2", RR_COLOR);
        ffOrientationColors.put("R2F1", RR_COLOR);
        //RR
        ffOrientationColors.put("R1F2", LL_COLOR);
        ffOrientationColors.put("F2R1", LL_COLOR);
        //RL
        ffOrientationColors.put("R1R2", RL_COLOR);
        ffOrientationColors.put("F2F1", RL_COLOR);
    }


    /**
     * Constructs ...
     */
    public AlignmentRenderer() {
        this.prefs = PreferenceManager.getInstance();

    }

    /**
     * Render a row of alignments in the given rectangle.
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
        Font font = FontManager.getFont(10);

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


                // If the alignment is 3 pixels or less,  draw alignment as a single block,
                // further detail would not be seen and just add to drawing overhead
                // Does the change for Bisulfite kill some machines?
                double pixelWidth = pixelEnd - pixelStart;
                if ((pixelWidth < 4) && !(AlignmentTrack.isBisulfiteColorType(renderOptions.getColorOption()) && (pixelWidth >= 1))) {
                    Color alignmentColor = getAlignmentColor(alignment, renderOptions);
                    Graphics2D g = context.getGraphic2DForColor(alignmentColor);
                    g.setFont(font);

                    int w = Math.max(1, (int) (pixelWidth));
                    int h = (int) Math.max(1, rect.getHeight() - 2);
                    int y = (int) (rect.getY() + (rect.getHeight() - h) / 2);
                    g.fillRect((int) pixelStart, y, w, h);
                } else if (alignment instanceof PairedAlignment) {
                    drawPairedAlignment((PairedAlignment) alignment, rect, context, renderOptions, leaveMargin, selectedReadNames, font);
                } else {
                    Color alignmentColor = getAlignmentColor(alignment, renderOptions);
                    Graphics2D g = context.getGraphic2DForColor(alignmentColor);
                    g.setFont(font);
                    drawAlignment(alignment, rect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames);
                }

            }

            // Optionally draw a border around the center base
            boolean showCenterLine = prefs.getAsBoolean(PreferenceManager.SAM_SHOW_CENTER_LINE);
            final int bottom = rect.y + rect.height;
            if (locScale < 5 && showCenterLine) {
                // Calculate center lines
                double center = (int) (context.getReferenceFrame().getCenter() - origin);
                int centerLeftP = (int) (center / locScale);
                int centerRightP = (int) ((center + 1) / locScale);
                //float transparency = Math.max(0.5f, (float) Math.round(10 * (1 - .75 * locScale)) / 10);
                Graphics2D gBlack = context.getGraphic2DForColor(Color.black); //new Color(0, 0, 0, transparency));
                GraphicUtils.drawDottedDashLine(gBlack, centerLeftP, rect.y, centerLeftP, bottom);
                if ((centerRightP - centerLeftP > 2)) {
                    GraphicUtils.drawDottedDashLine(gBlack, centerRightP, rect.y, centerRightP, bottom);
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

    /**
     * Draw a pair of alignments as a single "template".
     *
     * @param pair
     * @param rect
     * @param context
     * @param renderOptions
     * @param leaveMargin
     * @param selectedReadNames
     * @param font
     */
    private void drawPairedAlignment(
            PairedAlignment pair,
            Rectangle rect,
            RenderContext context,
            AlignmentTrack.RenderOptions renderOptions,
            boolean leaveMargin,
            Map<String, Color> selectedReadNames,
            Font font) {

        double locScale = context.getScale();
        Color alignmentColor = getAlignmentColor(pair.firstAlignment, renderOptions);
        Graphics2D g = context.getGraphic2DForColor(alignmentColor);
        g.setFont(font);
        drawAlignment(pair.firstAlignment, rect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames);

        if (pair.secondAlignment != null) {

            alignmentColor = getAlignmentColor(pair.secondAlignment, renderOptions);
            g = context.getGraphic2DForColor(alignmentColor);

            drawAlignment(pair.secondAlignment, rect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames);

            Graphics2D gLine = context.getGraphic2DForColor(grey1);
            double origin = context.getOrigin();
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
     * Draw a (possible) gapped alignment
     *
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

        // No blocks.  Note: SAM/BAM alignments always have at least 1 block
        if (blocks == null || blocks.length == 0) {
            drawSimpleAlignment(alignment, rect, g, context, renderOptions.flagUnmappedPairs);
            return;
        }

        // Get the terminal block (last block with respect to read direction).  This will have an "arrow" attached.
        AlignmentBlock terminalBlock = alignment.isNegativeStrand() ? blocks[0] : blocks[blocks.length - 1];


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


            // Get a graphics context for outlining reads
            Graphics2D outlineGraphics = context.getGraphic2DForColor(OUTLINE_COLOR);

            // Create polygon to represent the alignment.
            boolean isZeroQuality = alignment.getMappingQuality() == 0 && renderOptions.flagZeroQualityAlignments;

            // If we're zoomed in and this is a large block clip a pixel off each end.  TODO - why?
            if (highZoom && w > 10) {
                x++;
                w -= 2;
            }

            // If block is out of view skip -- this is important in the case of PacBio and other platforms with very long reads
            if (x + w >= rect.x && x <= rect.getMaxX()) {

                Shape blockShape = null;

                // If this is a terminal block draw the "arrow" to indicate strand position.  Otherwise draw a rectangle.
                if ((aBlock == terminalBlock) && w > 10 && h > 10) {

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
                    blockShape = new Polygon(xPoly, yPoly, xPoly.length);
                } else {
                    // Not a terminal block, or too small for arrow
                    blockShape = new Rectangle(x, y, w, h);
                }

                g.fill(blockShape);

                if (isZeroQuality) {
                    outlineGraphics.draw(blockShape);
                }

                if (renderOptions.flagUnmappedPairs && alignment.isPaired() && !alignment.getMate().isMapped()) {
                    Graphics2D cRed = context.getGraphic2DForColor(Color.red);
                    cRed.draw(blockShape);
                }

                if (selectedReadNames.containsKey(alignment.getReadName())) {
                    Color c = selectedReadNames.get(alignment.getReadName());
                    if (c == null) {
                        c = Color.blue;
                    }
                    Graphics2D cBlue = context.getGraphic2DForColor(c);
                    Stroke s = cBlue.getStroke();
                    cBlue.setStroke(thickStroke);
                    cBlue.draw(blockShape);
                    cBlue.setStroke(s);
                }

            }

            if ((locScale < 5) || (AlignmentTrack.isBisulfiteColorType(renderOptions.getColorOption()) && (locScale < 100))) // Is 100 here going to kill some machines? bpb
            {
                if (renderOptions.showMismatches || renderOptions.showAllBases) {
                    drawBases(context, rect, aBlock, alignmentColor, renderOptions);
                }
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
     * Draw bases for an alignment block.  The bases are "overlaid" on the block with a transparency value (alpha)
     * that is proportional to the base quality score.
     *
     * @param context
     * @param rect
     * @param block
     * @param alignmentColor
     */

    private void drawBases(RenderContext context,
                           Rectangle rect,
                           AlignmentBlock block,
                           Color alignmentColor,
                           RenderOptions renderOptions) {


        boolean shadeBases = renderOptions.shadeBases;
        ColorOption colorOption = renderOptions.getColorOption();

        // Disable showAllBases in bisulfite mode
        boolean showAllBases = renderOptions.showAllBases &&
                !(colorOption == ColorOption.BISULFITE || colorOption == ColorOption.NOMESEQ);

        double locScale = context.getScale();
        double origin = context.getOrigin();
        String chr = context.getChr();
        //String genomeId = context.getGenomeId();
        Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();

        byte[] read = block.getBases();
        boolean isSoftClipped = block.isSoftClipped();

        int start = block.getStart();
        int end = start + read.length;
        byte[] reference = isSoftClipped ? null : genome.getSequence(chr, start, end);


        if (read != null && read.length > 0 && reference != null) {

            // Compute bounds, get a graphics to use,  and compute a font
            int pY = (int) rect.getY();
            int dY = (int) rect.getHeight();
            int dX = (int) Math.max(1, (1.0 / locScale));
            Graphics2D g = (Graphics2D) context.getGraphics().create();
            if (dX >= 8) {
                Font f = FontManager.getFont(Font.BOLD, Math.min(dX, 12));
                g.setFont(f);
            }

            // Get the base qualities, start/end,  and reference sequence

            BisulfiteBaseInfo bisinfo = null;
            boolean nomeseqMode = (renderOptions.getColorOption().equals(AlignmentTrack.ColorOption.NOMESEQ));
            boolean bisulfiteMode = AlignmentTrack.isBisulfiteColorType(renderOptions.getColorOption());
            if (nomeseqMode) {
                bisinfo = new BisulfiteBaseInfoNOMeseq(reference, block, renderOptions.bisulfiteContext);
            } else if (bisulfiteMode) {
                bisinfo = new BisulfiteBaseInfo(reference, block, renderOptions.bisulfiteContext);
            }

            // Loop through base pair coordinates
            for (int loc = start; loc < end; loc++) {

                // Index into read array,  just the genomic location offset by
                // the start of this block
                int idx = loc - start;

                // Is this base a mismatch?  Note '=' means indicates a match by definition
                // If we do not have a valid reference we assume a match.  Soft clipped
                // bases are considered mismatched by definition
                boolean misMatch;
                if (isSoftClipped) {
                    misMatch = true;  // <= by definition, any matches are coincidence
                } else {
                    final byte refbase = reference[idx];
                    final byte readbase = read[idx];
                    misMatch = readbase != '=' &&
                            reference != null &&
                            idx < reference.length &&
                            refbase != 0 &&
                            !AlignmentUtils.compareBases(refbase, readbase);
                }


                if (showAllBases || (!bisulfiteMode && misMatch) ||
                        (bisulfiteMode && (!DisplayStatus.NOTHING.equals(bisinfo.getDisplayStatus(idx))))) {
                    char c = (char) read[idx];

                    Color color = nucleotideColors.get(c);
                    if (bisulfiteMode) color = bisinfo.getDisplayColor(idx);
                    if (color == null) {
                        color = Color.black;
                    }

                    if (shadeBases) {
                        byte qual = block.qualities[loc - start];
                        color = getShadedColor(qual, color, alignmentColor, prefs);
                    }

                    double bisulfiteXaxisShift = (bisulfiteMode) ? bisinfo.getXaxisShift(idx) : 0;

                    // If there is room for text draw the character, otherwise
                    // just draw a rectangle to represent the
                    int pX0 = (int) (((double) loc + bisulfiteXaxisShift - (double) origin) / (double) locScale);

                    // Don't draw out of clipping rect
                    if (pX0 > rect.getMaxX()) {
                        break;
                    } else if (pX0 + dX < rect.getX()) {
                        continue;
                    }


                    BisulfiteBaseInfo.DisplayStatus bisstatus = (bisinfo == null) ? null : bisinfo.getDisplayStatus(idx);
                    // System.err.printf("Draw text?  dY=%d, dX=%d, bismode=%s, dispStatus=%s\n",dY,dX,!bisulfiteMode || bisulfiteMode,bisstatus);
                    if (((dY >= 12) && (dX >= 8)) && (!bisulfiteMode || (bisulfiteMode && bisstatus.equals(DisplayStatus.CHARACTER)))) {
                        g.setColor(color);
                        GraphicUtils.drawCenteredText(g, new char[]{c}, pX0, pY + 1, dX, dY - 2);
                    } else {

                        int pX0i = pX0, dXi = dX;

                        // If bisulfite mode, we expand the rectangle to make it more visible
                        if (bisulfiteMode && bisstatus.equals(DisplayStatus.COLOR)) {
                            if (dXi < 3) {
                                int expansion = dXi;
                                pX0i -= expansion;
                                dXi += (2 * expansion);
                            }
                        }

                        int dW = (dXi > 4 ? dXi - 1 : dXi);

                        if (color != null) {
                            g.setColor(color);
                            if (dY < 10) {
                                g.fillRect(pX0i, pY, dXi, dY);
                            } else {
                                g.fillRect(pX0i, pY + 1, dW, dY - 3);
                            }
                        }
                    }
                }

            }
        }
    }

    private Color getShadedColor(byte qual, Color foregroundColor, Color backgroundColor, PreferenceManager prefs) {
        float alpha = 0;
        int minQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN);
        if (qual < minQ) {
            alpha = 0.1f;
        } else {
            int maxQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX);
            alpha = Math.max(0.1f, Math.min(1.0f, 0.1f + 0.9f * (qual - minQ) / (maxQ - minQ)));
        }
        // Round alpha to nearest 0.1
        alpha = ((int) (alpha * 10 + 0.5f)) / 10.0f;

        if (alpha >= 1) {
            return foregroundColor;
        }
        Color color = ColorUtilities.getCompositeColor(backgroundColor, foregroundColor, alpha);
        return color;
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

    private Color getAlignmentColor(Alignment alignment, AlignmentTrack.RenderOptions renderOptions) {

        // Set color used to draw the feature.  Highlight features that intersect the
        // center line.  Also restorePersistentState row "score" if alignment intersects center line

        String lb = alignment.getLibrary();
        if (lb == null) lb = "null";
        PEStats peStats = renderOptions.peStats.get(lb);

        Color c = alignment.getDefaultColor();
        switch (renderOptions.getColorOption()) {
            case BISULFITE:
                // Just a simple forward/reverse strand color scheme that won't clash with the
                // methylation rectangles.
                c = (alignment.getFirstOfPairStrand() == Strand.POSITIVE) ? bisulfiteColorFw1 : bisulfiteColorRev1;

//                if (alignment.isNegativeStrand()) {
//                    c = (alignment.isSecondOfPair()) ? bisulfiteColorRev2 : bisulfiteColorRev1;
//                } else {
//                    c = (alignment.isSecondOfPair()) ? bisulfiteColorFw2 : bisulfiteColorFw1;
//                }
                break;
            case NOMESEQ:
                c = nomeseqColor;
                break;

            case INSERT_SIZE:
                boolean isPairedAlignment = alignment instanceof PairedAlignment;
                if (alignment.isPaired() && alignment.getMate().isMapped() || isPairedAlignment) {
                    boolean sameChr = isPairedAlignment ||
                            alignment.getMate().getChr().equals(alignment.getChr());
                    if (sameChr) {
                        int readDistance = Math.abs(alignment.getInferredInsertSize());
                        if (readDistance > 0) {
                            int minThreshold = renderOptions.getMinInsertSize();
                            int maxThreshold = renderOptions.getMaxInsertSize();
                            if (renderOptions.isComputeIsizes() && renderOptions.peStats != null) {
                                if (peStats != null) {
                                    minThreshold = peStats.getMinThreshold();
                                    maxThreshold = peStats.getMaxThreshold();
                                }

                            }

                            if (readDistance < minThreshold) {
                                c = smallISizeColor;
                            } else if (readDistance > maxThreshold) {
                                c = largeISizeColor;
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
                c = getOrientationColor(alignment, peStats);
                break;
            case READ_STRAND:
                if (alignment.isNegativeStrand()) {
                    c = negStrandColor;
                } else {
                    c = posStrandColor;
                }
                break;
            case FRAGMENT_STRAND:
                final Strand fragmentStrand = alignment.getFirstOfPairStrand();
                if (fragmentStrand == Strand.NEGATIVE) {
                    c = negStrandColor;
                } else if (fragmentStrand == Strand.POSITIVE) {
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
            case TAG:
                final String tag = renderOptions.getColorByTag();
                if (tag != null) {
                    Object tagValue = alignment.getAttribute(tag);
                    if (tagValue != null) {
                        c = readGroupColors.get(tagValue);
                        if (c == null) {
                            c = ColorUtilities.randomColor(readGroupColors.size() + 1);
                            readGroupColors.put(tagValue.toString(), c);
                        }
                    }
                }
                break;

            default:
//                if (renderOptions.shadeCenters && center >= alignment.getStart() && center <= alignment.getEnd()) {
//                    if (locScale < 1) {
//                        c = grey2;
//                    }
//                }

        }
        if (c == null) c = grey1;

        if (alignment.getMappingQuality() == 0 && renderOptions.flagZeroQualityAlignments) {
            // Maping Q = 0
            float alpha = 0.15f;
            // Assuming white background TODO -- this should probably be passed in
            return ColorUtilities.getCompositeColor(Color.white, c, alpha);
        }

        return c;
    }


    /**
     * @return
     */

    private Color getOrientationColor(Alignment alignment, PEStats peStats) {

        Color c = null;
        if (alignment.isPaired() && !alignment.isProperPair()) {

            final String pairOrientation = alignment.getPairOrientation();
            if (peStats != null) {
                PEStats.Orientation libraryOrientation = peStats.getOrientation();
                switch (libraryOrientation) {
                    case FR:
                        if (!alignment.isSmallInsert()) {
                            // if the isize < read length the reads might overlap, invalidating this test
                            c = frOrientationColors.get(pairOrientation);
                        }
                        break;
                    case RF:
                        c = rfOrientationColors.get(pairOrientation);
                        break;
                    case FF:
                        c = ffOrientationColors.get(pairOrientation);
                        break;
                }

            } else {
                // No peStats for this library
                if (alignment.getAttribute("CS") != null) {
                    c = ffOrientationColors.get(pairOrientation);
                } else {
                    c = frOrientationColors.get(pairOrientation);
                }
            }
        }

        return c == null ? grey1 : c;

    }

    /**
     * @return
     */


    public static Map<Character, Color> getNucleotideColors() {
        return nucleotideColors;
    }
}
