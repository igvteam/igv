/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.sam.AlignmentTrack.ColorOption;
import org.broad.igv.sam.AlignmentTrack.RenderOptions;
import org.broad.igv.sam.AlignmentTrack.ShadeBasesOption;
import org.broad.igv.sam.BisulfiteBaseInfo.DisplayStatus;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.color.ColorPalette;
import org.broad.igv.ui.color.ColorTable;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.util.ChromosomeColors;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.QuadCurve2D;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AlignmentRenderer implements FeatureRenderer {

    private static Logger log = Logger.getLogger(AlignmentRenderer.class);

    public static final Color GROUP_DIVIDER_COLOR = new Color(200, 200, 200);

    // A "dummy" reference for soft-clipped reads.
    private static byte[] softClippedReference = new byte[1000];

    private static Color smallISizeColor = new Color(0, 0, 150);
    private static Color largeISizeColor = new Color(150, 0, 0);
    private static Color purple = new Color(118, 24, 220);
    private static Color deletionColor = Color.black;
    private static Color skippedColor = new Color(150, 184, 200);
    public static Color grey1 = new Color(200, 200, 200);

    private static Stroke thickStroke = new BasicStroke(2.0f);

    // Bisulfite constants
    private final Color bisulfiteColorFw1 = new Color(195, 195, 195);
    private final Color bisulfiteColorRev1 = new Color(195, 210, 195);
    private final Color nomeseqColor = new Color(195, 195, 195);

    public static final Color negStrandColor = new Color(150, 150, 230);
    public static final Color posStrandColor = new Color(230, 150, 150);

    private ColorTable readGroupColors;
    private ColorTable sampleColors;
    private ColorTable tagValueColors;

    private final Color LR_COLOR = grey1; // "Normal" alignment color
    //private final Color LR_COLOR_12 = new Color(190, 190, 210);
    //private final Color LR_COLOR_21 = new Color(210, 190, 190);
    private final Color RL_COLOR = new Color(0, 150, 0);
    private final Color RR_COLOR = new Color(20, 50, 200);
    private final Color LL_COLOR = new Color(0, 150, 150);

    private final Color OUTLINE_COLOR = new Color(185, 185, 185);

    private static Map<String, AlignmentTrack.OrientationType> frOrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> f1f2OrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> f2f1OrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> rfOrientationTypes;
    private Map<AlignmentTrack.OrientationType, Color> typeToColorMap;

    static {
        initializeTagTypes();
    }

    PreferenceManager prefs;

    private static AlignmentRenderer instance;

    private TreeSet<Shape> arcsByStart;
    private TreeSet<Shape> arcsByEnd;
    private HashMap<Shape, Alignment> curveMap;

    public static AlignmentRenderer getInstance() {
        if (instance == null) {
            instance = new AlignmentRenderer();
        }
        return instance;
    }


    private AlignmentRenderer() {
        this.prefs = PreferenceManager.getInstance();
        initializeTagColors();
        curveMap = new HashMap<Shape, Alignment>();

        arcsByStart = new TreeSet<Shape>(new Comparator<Shape>() {

            public int compare(Shape o1, Shape o2) {
                double x1 = o1.getBounds().getMinX();
                double x2 = o2.getBounds().getMinX();
                return (int) Math.signum(x1 - x2);
            }
        });

        arcsByEnd = new TreeSet<Shape>(new Comparator<Shape>() {

            public int compare(Shape o1, Shape o2) {
                double x1 = o1.getBounds().getMaxX();
                double x2 = o2.getBounds().getMaxX();
                return (int) Math.signum(x1 - x2);
            }
        });
    }

    private static void initializeTagTypes() {
        // pre-seed from orientation colors

        // fr Orientations (e.g. Illumina paired-end libraries)
        frOrientationTypes = new HashMap();
        //LR
        frOrientationTypes.put("F1R2", AlignmentTrack.OrientationType.LR);
        frOrientationTypes.put("F2R1", AlignmentTrack.OrientationType.LR);
        frOrientationTypes.put("F R ", AlignmentTrack.OrientationType.LR);
        frOrientationTypes.put("FR", AlignmentTrack.OrientationType.LR);
        //LL
        frOrientationTypes.put("F1F2", AlignmentTrack.OrientationType.LL);
        frOrientationTypes.put("F2F1", AlignmentTrack.OrientationType.LL);
        frOrientationTypes.put("F F ", AlignmentTrack.OrientationType.LL);
        frOrientationTypes.put("FF", AlignmentTrack.OrientationType.LL);
        //RR
        frOrientationTypes.put("R1R2", AlignmentTrack.OrientationType.RR);
        frOrientationTypes.put("R2R1", AlignmentTrack.OrientationType.RR);
        frOrientationTypes.put("R R ", AlignmentTrack.OrientationType.RR);
        frOrientationTypes.put("RR", AlignmentTrack.OrientationType.RR);
        //RL
        frOrientationTypes.put("R1F2", AlignmentTrack.OrientationType.RL);
        frOrientationTypes.put("R2F1", AlignmentTrack.OrientationType.RL);
        frOrientationTypes.put("R F ", AlignmentTrack.OrientationType.RL);
        frOrientationTypes.put("RF", AlignmentTrack.OrientationType.RL);

        // rf orienation  (e.g. Illumina mate-pair libraries)
        rfOrientationTypes = new HashMap();
        //LR
        rfOrientationTypes.put("R1F2", AlignmentTrack.OrientationType.LR);
        rfOrientationTypes.put("R2F1", AlignmentTrack.OrientationType.LR);
        //rfOrientationTypes.put("R F ", AlignmentTrack.OrientationType.LR);
        //rfOrientationTypes.put("RF", AlignmentTrack.OrientationType.LR);
        //LL
        rfOrientationTypes.put("R1R2", AlignmentTrack.OrientationType.LL);
        rfOrientationTypes.put("R2R1", AlignmentTrack.OrientationType.LL);
        rfOrientationTypes.put("R R ", AlignmentTrack.OrientationType.LL);
        rfOrientationTypes.put("RR ", AlignmentTrack.OrientationType.LL);

        rfOrientationTypes.put("F1F2", AlignmentTrack.OrientationType.RR);
        rfOrientationTypes.put("F2F1", AlignmentTrack.OrientationType.RR);
        rfOrientationTypes.put("F F ", AlignmentTrack.OrientationType.RR);
        rfOrientationTypes.put("FF", AlignmentTrack.OrientationType.RR);
        //RL
        rfOrientationTypes.put("F1R2", AlignmentTrack.OrientationType.RL);
        rfOrientationTypes.put("F2R1", AlignmentTrack.OrientationType.RL);
        rfOrientationTypes.put("F R ", AlignmentTrack.OrientationType.RL);
        rfOrientationTypes.put("FR", AlignmentTrack.OrientationType.RL);

        // f1f2 orienation  (e.g. SOLID libraries, AlignmentTrack.OrientationType.second read appears first on + strand (leftmost))
        f2f1OrientationTypes = new HashMap();
        //LR
        f2f1OrientationTypes.put("F2F1", AlignmentTrack.OrientationType.LR);
        f2f1OrientationTypes.put("R1R2", AlignmentTrack.OrientationType.LR);

        //LL
        f2f1OrientationTypes.put("F2R1", AlignmentTrack.OrientationType.LL);
        f2f1OrientationTypes.put("R1F2", AlignmentTrack.OrientationType.LL);

        //RR
        f2f1OrientationTypes.put("R2F1", AlignmentTrack.OrientationType.RR);
        f2f1OrientationTypes.put("F1R2", AlignmentTrack.OrientationType.RR);

        //RL
        f2f1OrientationTypes.put("R2R1", AlignmentTrack.OrientationType.RL);
        f2f1OrientationTypes.put("F1F2", AlignmentTrack.OrientationType.RL);

        // f1f2 orienation  (e.g. SOLID libraries, AlignmentTrack.OrientationType.actually is this one even possible?)
        f1f2OrientationTypes = new HashMap();
        //LR
        f1f2OrientationTypes.put("F1F2", AlignmentTrack.OrientationType.LR);
        f1f2OrientationTypes.put("R2R1", AlignmentTrack.OrientationType.LR);
        //LL
        f1f2OrientationTypes.put("F1R2", AlignmentTrack.OrientationType.LL);
        f1f2OrientationTypes.put("R2F1", AlignmentTrack.OrientationType.LL);
        //RR
        f1f2OrientationTypes.put("R1F2", AlignmentTrack.OrientationType.RR);
        f1f2OrientationTypes.put("F2R1", AlignmentTrack.OrientationType.RR);
        //RL
        f1f2OrientationTypes.put("R1R2", AlignmentTrack.OrientationType.RL);
        f1f2OrientationTypes.put("F2F1", AlignmentTrack.OrientationType.RL);
    }

    private void initializeTagColors() {
        ColorPalette palette = ColorUtilities.getPalette("Pastel 1");  // TODO let user choose
        readGroupColors = new PaletteColorTable(palette);
        sampleColors = new PaletteColorTable(palette);
        tagValueColors = new PaletteColorTable(palette);

        typeToColorMap = new HashMap<AlignmentTrack.OrientationType, Color>(5);
        typeToColorMap.put(AlignmentTrack.OrientationType.LL, LL_COLOR);
        typeToColorMap.put(AlignmentTrack.OrientationType.LR, LR_COLOR);
        typeToColorMap.put(AlignmentTrack.OrientationType.RL, RL_COLOR);
        typeToColorMap.put(AlignmentTrack.OrientationType.RR, RR_COLOR);
        typeToColorMap.put(null, grey1);
    }

    /**
     * Render a row of alignments in the given rectangle.
     */
    public void renderAlignments(List<Alignment> alignments,
                                 RenderContext context,
                                 Rectangle rowRect,
                                 Rectangle trackRect, RenderOptions renderOptions,
                                 boolean leaveMargin,
                                 Map<String, Color> selectedReadNames) {

        double origin = context.getOrigin();
        double locScale = context.getScale();
        Font font = FontManager.getFont(10);

        if ((alignments != null) && (alignments.size() > 0)) {

            int lastPixelDrawn = -1;

            for (Alignment alignment : alignments) {
                // Compute the start and dend of the alignment in pixels
                double pixelStart = ((alignment.getStart() - origin) / locScale);
                double pixelEnd = ((alignment.getEnd() - origin) / locScale);

                // If the any part of the feature fits in the track rectangle draw  it
                if (pixelEnd < rowRect.x) {
                    continue;
                } else if (pixelStart > rowRect.getMaxX()) {
                    break;
                }


                // If the alignment is 3 pixels or less,  draw alignment as a single block,
                // further detail would not be seen and just add to drawing overhead
                // Does the change for Bisulfite kill some machines?
                double pixelWidth = pixelEnd - pixelStart;
                if ((pixelWidth < 4) && !(AlignmentTrack.isBisulfiteColorType(renderOptions.getColorOption()) && (pixelWidth >= 1))) {

                    Color alignmentColor = getAlignmentColor(alignment, renderOptions);

                    // Optimization for really zoomed out views.  If this alignment occupies screen space already taken,
                    // and it is the default color, skip drawing.
                    if (pixelEnd <= lastPixelDrawn && alignmentColor == alignment.getDefaultColor()) {
                        continue;
                    }


                    Graphics2D g = context.getGraphic2DForColor(alignmentColor);
                    g.setFont(font);

                    int w = Math.max(1, (int) (pixelWidth));
                    int h = (int) Math.max(1, rowRect.getHeight() - 2);
                    int y = (int) (rowRect.getY() + (rowRect.getHeight() - h) / 2);
                    g.fillRect((int) pixelStart, y, w, h);
                    lastPixelDrawn = (int) pixelStart + w;
                } else if (alignment instanceof PairedAlignment) {
                    drawPairedAlignment((PairedAlignment) alignment, rowRect, trackRect, context, renderOptions, leaveMargin, selectedReadNames, font);
                } else {
                    Color alignmentColor = getAlignmentColor(alignment, renderOptions);
                    Graphics2D g = context.getGraphic2DForColor(alignmentColor);
                    g.setFont(font);
                    drawAlignment(alignment, rowRect, trackRect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames);
                }
            }

            // Optionally draw a border around the center base
            boolean showCenterLine = prefs.getAsBoolean(PreferenceManager.SAM_SHOW_CENTER_LINE);
            final int bottom = rowRect.y + rowRect.height;
            if (locScale < 5 && showCenterLine) {
                // Calculate center lines
                double center = (int) (context.getReferenceFrame().getCenter() - origin);
                int centerLeftP = (int) (center / locScale);
                int centerRightP = (int) ((center + 1) / locScale);
                //float transparency = Math.max(0.5f, (float) Math.round(10 * (1 - .75 * locScale)) / 10);
                Graphics2D gBlack = context.getGraphic2DForColor(Color.black); //new Color(0, 0, 0, transparency));
                GraphicUtils.drawDottedDashLine(gBlack, centerLeftP, rowRect.y, centerLeftP, bottom);
                if ((centerRightP - centerLeftP > 2)) {
                    GraphicUtils.drawDottedDashLine(gBlack, centerRightP, rowRect.y, centerRightP, bottom);
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
     * @param rowRect
     * @param context
     * @param renderOptions
     * @param leaveMargin
     * @param selectedReadNames
     * @param font
     */
    private void drawPairedAlignment(
            PairedAlignment pair,
            Rectangle rowRect,
            Rectangle trackRect,
            RenderContext context,
            AlignmentTrack.RenderOptions renderOptions,
            boolean leaveMargin,
            Map<String, Color> selectedReadNames,
            Font font) {

        //Only plot outliers
        if (renderOptions.isPairedArcView() && getOutlierStatus(pair, renderOptions) == 0) {
            return;
        }

        double locScale = context.getScale();

        Color alignmentColor1;
        Color alignmentColor2 = null;
        if (renderOptions.isPairedArcView()) {
            renderOptions.setColorOption(ColorOption.INSERT_SIZE);
            alignmentColor1 = getAlignmentColor(pair, renderOptions);
            alignmentColor2 = alignmentColor1;
        } else {
            alignmentColor1 = getAlignmentColor(pair.firstAlignment, renderOptions);
        }

        Graphics2D g = context.getGraphic2DForColor(alignmentColor1);
        g.setFont(font);
        drawAlignment(pair.firstAlignment, rowRect, trackRect, g, context, alignmentColor1, renderOptions, leaveMargin, selectedReadNames);

        //If the paired alignment is in memory, we draw it.
        //However, we get the coordinates from the first alignment
        if (pair.secondAlignment != null) {

            if (alignmentColor2 == null) {
                alignmentColor2 = getAlignmentColor(pair.secondAlignment, renderOptions);
            }
            g = context.getGraphic2DForColor(alignmentColor2);

            drawAlignment(pair.secondAlignment, rowRect, trackRect, g, context, alignmentColor2, renderOptions, leaveMargin, selectedReadNames);
        } else {
            return;
        }

        Color lineColor = grey1;

        if (alignmentColor1.equals(alignmentColor2) || pair.secondAlignment == null) {
            lineColor = alignmentColor1;
        }
        Graphics2D gLine = context.getGraphic2DForColor(lineColor);

        double origin = context.getOrigin();
        int startX = (int) ((pair.firstAlignment.getEnd() - origin) / locScale);
        int endX = (int) ((pair.firstAlignment.getMate().getStart() - origin) / locScale);

        int h = (int) Math.max(1, rowRect.getHeight() - (leaveMargin ? 2 : 0));
        int y = (int) (rowRect.getY());


        if (renderOptions.isPairedArcView()) {
            int relation = compareToBounds(pair, renderOptions);
            if (relation <= -1 || relation >= +1) {
                return;
            }
            GeneralPath path = new GeneralPath(GeneralPath.WIND_NON_ZERO, 4);
            int curveHeight = (int) Math.log(endX - startX) * h;

            double botY = y + h / 2;
            double topY = y + h / 2 - curveHeight;
            double midX = (endX + startX) / 2;

            path.moveTo(startX, botY);
            path.quadTo(midX, topY, endX, botY);
            path.quadTo(midX, topY - 2, startX, botY);
            path.closePath();
            arcsByStart.add(path);
            arcsByEnd.add(path);
            curveMap.put(path, pair);
            gLine.setColor(alignmentColor2);

            gLine.draw(path);
        } else {
            startX = Math.max(rowRect.x, startX);
            endX = Math.min(rowRect.x + rowRect.width, endX);
            gLine.drawLine(startX, y + h / 2, endX, y + h / 2);
        }

    }

    /**
     * Draw a (possibly gapped) alignment
     *
     * @param alignment
     * @param rowRect
     * @param trackRect
     * @param g
     * @param context
     * @param alignmentColor
     * @param renderOptions
     * @param leaveMargin
     * @param selectedReadNames
     */
    private void drawAlignment(
            Alignment alignment,
            Rectangle rowRect,
            Rectangle trackRect,
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
            drawSimpleAlignment(alignment, rowRect, g, context, renderOptions.flagUnmappedPairs);
            return;
        }


        // Get the terminal block (last block with respect to read direction).  This will have an "arrow" attached.
        AlignmentBlock terminalBlock = alignment.isNegativeStrand() ? blocks[0] : blocks[blocks.length - 1];

        int lastBlockEnd = Integer.MIN_VALUE;

        int blockNumber = -1;
        char[] gapTypes = alignment.getGapTypes();
        boolean highZoom = locScale < 0.1251;

        // Get a graphics context for outlining reads
        Graphics2D outlineGraphics = context.getGraphic2DForColor(OUTLINE_COLOR);
        Graphics2D terminalGrpahics = context.getGraphic2DForColor(Color.DARK_GRAY);

        boolean isZeroQuality = alignment.getMappingQuality() == 0 && renderOptions.flagZeroQualityAlignments;
        int h = (int) Math.max(1, rowRect.getHeight() - (leaveMargin ? 2 : 0));
        int y = (int) (rowRect.getY());

        for (AlignmentBlock aBlock : alignment.getAlignmentBlocks()) {
            blockNumber++;
            int blockPixelStart = (int) ((aBlock.getStart() - origin) / locScale);
            int blockPixelWidth = (int) Math.ceil(aBlock.getBases().length / locScale);

            // If we're zoomed in and this is a large block clip a pixel off each end.  TODO - why?
            if (highZoom && blockPixelWidth > 10) {
                blockPixelStart++;
                blockPixelWidth -= 2;
            }

            // If block is out of view skip -- this is important in the case of PacBio and other platforms with very long reads
            if (blockPixelStart + blockPixelWidth >= rowRect.x && blockPixelStart <= rowRect.getMaxX()) {

                Shape blockShape = null;

                // If this is a terminal block draw the "arrow" to indicate strand position.  Otherwise draw a rectangle.
                if ((aBlock == terminalBlock) && blockPixelWidth > 10)
                    if (h > 10) {

                        int arrowLength = Math.min(5, blockPixelWidth / 6);

                        // Don't draw off edge of clipping rect
                        if (blockPixelStart < rowRect.x && (blockPixelStart + blockPixelWidth) > (rowRect.x + rowRect.width)) {
                            blockPixelStart = rowRect.x;
                            blockPixelWidth = rowRect.width;
                            arrowLength = 0;
                        } else if (blockPixelStart < rowRect.x) {
                            int delta = rowRect.x - blockPixelStart;
                            blockPixelStart = rowRect.x;
                            blockPixelWidth -= delta;
                            if (alignment.isNegativeStrand()) {
                                arrowLength = 0;
                            }
                        } else if ((blockPixelStart + blockPixelWidth) > (rowRect.x + rowRect.width)) {
                            blockPixelWidth -= ((blockPixelStart + blockPixelWidth) - (rowRect.x + rowRect.width));
                            if (!alignment.isNegativeStrand()) {
                                arrowLength = 0;
                            }
                        }

                        int[] xPoly;
                        int[] yPoly = {y, y, y + h / 2, y + h, y + h};

                        if (alignment.isNegativeStrand()) {
                            xPoly = new int[]{blockPixelStart + blockPixelWidth, blockPixelStart, blockPixelStart - arrowLength, blockPixelStart, blockPixelStart + blockPixelWidth};
                        } else {
                            xPoly = new int[]{blockPixelStart, blockPixelStart + blockPixelWidth, blockPixelStart + blockPixelWidth + arrowLength, blockPixelStart + blockPixelWidth, blockPixelStart};
                        }
                        blockShape = new Polygon(xPoly, yPoly, xPoly.length);
                    } else {
                        // Terminal block, but not enough height for arrow.  Indicate with a line
                        int tH = Math.max(1, h - 1);
                        if (alignment.isNegativeStrand()) {
                            blockShape = new Rectangle(blockPixelStart, y, blockPixelWidth, h);
                            terminalGrpahics.drawLine(blockPixelStart, y, blockPixelStart, y + tH);
                        } else {
                            blockShape = new Rectangle(blockPixelStart, y, blockPixelWidth, h);
                            terminalGrpahics.drawLine(blockPixelStart + blockPixelWidth + 1, y, blockPixelStart + blockPixelWidth + 1, y + tH);
                        }
                    }
                else {
                    // Not a terminal block, or too small for arrow
                    blockShape = new Rectangle(blockPixelStart, y, blockPixelWidth, h);
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
                    drawBases(context, rowRect, aBlock, alignmentColor, renderOptions);
                }
            }

            // Draw connecting lines between blocks, if in view
            if (lastBlockEnd > Integer.MIN_VALUE && blockPixelStart > rowRect.x) {
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

                int startX = Math.max(rowRect.x, lastBlockEnd);
                int endX = Math.min(rowRect.x + rowRect.width, blockPixelStart);

                gLine.drawLine(startX, y + h / 2, endX, y + h / 2);
                gLine.setStroke(stroke);
            }
            lastBlockEnd = blockPixelStart + blockPixelWidth;

            // Next block cannot start before lastBlockEnd.  If its out of view we are done.
            if (lastBlockEnd > rowRect.getMaxX()) {
                break;
            }

        }

        // Render insertions if locScale ~ 0.25 (base level)
        if (locScale < 0.25) {
            drawInsertions(origin, rowRect, locScale, alignment, context);
        }


        //Draw straight line up for viewing arc pairs, if mate on a different chromosome
        if (renderOptions.isPairedArcView()) {
            try {
                Graphics2D gLine = context.getGraphic2DForColor(alignmentColor);
                if (!alignment.getChr().equalsIgnoreCase(alignment.getMate().getChr())) {
                    gLine.drawLine(lastBlockEnd, y + h / 2, lastBlockEnd, (int) trackRect.getMinY());
                }

            } catch (NullPointerException e) {
                //Don't have the info, don't plot anything
            }
        }

    }


    /**
     * Draw bases for an alignment block.  The bases are "overlaid" on the block with a transparency value (alpha)
     * that is proportional to the base quality score, or flow signal deviation, whichever is selected.
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


        ShadeBasesOption shadeBasesOption = renderOptions.shadeBasesOption;
        ColorOption colorOption = renderOptions.getColorOption();

        // Disable showAllBases in bisulfite mode
        boolean showAllBases = renderOptions.showAllBases &&
                !(colorOption == ColorOption.BISULFITE || colorOption == ColorOption.NOMESEQ);

        double locScale = context.getScale();
        double origin = context.getOrigin();
        String chr = context.getChr();
        //String genomeId = context.getGenomeId();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        byte[] read = block.getBases();
        boolean isSoftClipped = block.isSoftClipped();

        final int start = block.getStart();
        final int end = start + read.length;
        final byte[] reference = isSoftClipped ? softClippedReference : genome.getSequence(chr, start, end);

        if (read != null && read.length > 0 && reference != null) {

            // Compute bounds, get a graphics to use,  and compute a font
            int pY = (int) rect.getY();
            int dY = (int) rect.getHeight();
            int dX = (int) Math.max(1, (1.0 / locScale));
            Graphics2D g = (Graphics2D) context.getGraphics().create();
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
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
                // If we do not have a valid reference we assume a match.
                boolean misMatch = false;
                if (isSoftClipped) {
                    // Goby will return '=' characters when the soft-clip happens to match the reference.
                    // It could actually be useful to see which part of the soft clipped bases match, to help detect
                    // cases when an aligner clipped too much.
                    final byte readbase = read[idx];
                    misMatch = readbase != '=';  // mismatch, except when the soft-clip has an '=' base.
                } else {
                    if (reference != null) {
                        final int referenceLength = reference.length;
                        final byte refbase = idx < referenceLength ? reference[idx] : 0;
                        final byte readbase = read[idx];
                        misMatch = readbase != '=' &&
                                reference != null &&
                                idx < referenceLength &&
                                refbase != 0 &&
                                !AlignmentUtils.compareBases(refbase, readbase);
                    }
                }

                if (showAllBases || (!bisulfiteMode && misMatch) ||
                        (bisulfiteMode && (!DisplayStatus.NOTHING.equals(bisinfo.getDisplayStatus(idx))))) {
                    char c = (char) read[idx];

                    Color color = Globals.nucleotideColors.get(c);
                    if (bisulfiteMode) color = bisinfo.getDisplayColor(idx);
                    if (color == null) {
                        color = Color.black;
                    }

                    if (ShadeBasesOption.QUALITY == shadeBasesOption) {
                        byte qual = block.qualities[loc - start];
                        color = getShadedColor(qual, color, alignmentColor, prefs);
                    } else if (ShadeBasesOption.FLOW_SIGNAL_DEVIATION_READ == shadeBasesOption || ShadeBasesOption.FLOW_SIGNAL_DEVIATION_REFERENCE == shadeBasesOption) {
                        if (block.hasFlowSignals()) {
                            int flowSignal = (int) block.getFlowSignalSubContext(loc - start).getCurrentSignal();
                            int expectedFlowSignal;
                            if (ShadeBasesOption.FLOW_SIGNAL_DEVIATION_READ == shadeBasesOption) {
                                expectedFlowSignal = 100 * (short) ((flowSignal + 50.0) / 100.0);
                            } else {
                                // NB: this may estimate the reference homopolymer length incorrect in some cases, especially when we have
                                // an overcall/undercall situation.  Proper estimation of the reads observed versus expected homopolymer
                                // length should use flow signal alignment (SamToFlowspace): https://github.com/iontorrent/Ion-Variant-Hunter
                                if (!misMatch) {
                                    final byte readbase = read[idx];
                                    byte refbase = reference[idx];
                                    int pos; // zero based
                                    expectedFlowSignal = 100;

                                    // Count HP length
                                    pos = start + idx - 1;
                                    while (0 <= pos && genome.getReference(chr, pos) == refbase) {
                                        pos--;
                                        expectedFlowSignal += 100;
                                    }
                                    pos = start + idx + 1;
                                    while (pos < genome.getChromosome(chr).getLength() && genome.getReference(chr, pos) == refbase) {
                                        pos++;
                                        expectedFlowSignal += 100;
                                    }
                                } else {
                                    expectedFlowSignal = 0;
                                }
                            }
                            int flowSignalDiff = (expectedFlowSignal < flowSignal) ? (flowSignal - expectedFlowSignal) : (expectedFlowSignal - flowSignal);
                            // NB: this next section is some mangling to use the base quality color preferences...
                            if (flowSignalDiff <= 0) {
                                flowSignalDiff = 0;
                            } else if (50 < flowSignalDiff) {
                                flowSignalDiff = 50;
                            }
                            flowSignalDiff = 50 - flowSignalDiff; // higher is better
                            int minQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN);
                            int maxQ = prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX);
                            flowSignalDiff = flowSignalDiff * (maxQ - minQ) / 50;
                            byte qual;
                            int pos = start + idx;
                            if (Byte.MAX_VALUE < flowSignalDiff) {
                                qual = Byte.MAX_VALUE;
                            } else if (flowSignalDiff < Byte.MIN_VALUE) {
                                qual = Byte.MIN_VALUE;
                            } else {
                                qual = (byte) flowSignalDiff;
                            }
                            // Finally, get the color
                            color = getShadedColor(qual, color, alignmentColor, prefs);
                        }
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


        Color defaultColor = alignment.getDefaultColor();
        Color c = defaultColor;
        ColorOption colorOption = renderOptions.getColorOption();
        switch (colorOption) {
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

            case UNEXPECTED_PAIR:
            case PAIR_ORIENTATION:
                c = getOrientationColor(alignment, getPEStats(alignment, renderOptions));
                if (colorOption == ColorOption.PAIR_ORIENTATION || c != defaultColor) {
                    break;
                }
            case INSERT_SIZE:
                boolean isPairedAlignment = alignment instanceof PairedAlignment;
                if ((alignment.isPaired() && alignment.getMate().isMapped()) || isPairedAlignment) {
                    boolean sameChr = isPairedAlignment || alignment.getMate().getChr().equals(alignment.getChr());
                    if (sameChr) {
                        int readDistance = Math.abs(alignment.getInferredInsertSize());
                        if (readDistance != 0) {

                            int minThreshold = renderOptions.getMinInsertSize();
                            int maxThreshold = renderOptions.getMaxInsertSize();
                            PEStats peStats = getPEStats(alignment, renderOptions);
                            if (renderOptions.isComputeIsizes() && peStats != null) {
                                minThreshold = peStats.getMinThreshold();
                                maxThreshold = peStats.getMaxThreshold();
                            }

                            if (readDistance < minThreshold) {
                                c = smallISizeColor;
                            } else if (readDistance > maxThreshold) {
                                c = largeISizeColor;
                            }
                        }
                        //return renderOptions.insertSizeColorScale.getColor(readDistance);
                    } else {
                        c = ChromosomeColors.getColor(alignment.getMate().getChr());
                        if (c == null) {
                            c = Color.black;
                        }
                    }
                }


                break;
            case READ_STRAND:
                if (alignment.isNegativeStrand()) {
                    c = negStrandColor;
                } else {
                    c = posStrandColor;
                }
                break;
            case FIRST_OF_PAIR_STRAND:
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
                }
                break;
            case SAMPLE:
                String sample = alignment.getSample();
                if (sample != null) {
                    c = sampleColors.get(sample);
                }
                break;
            case TAG:
                final String tag = renderOptions.getColorByTag();
                if (tag != null) {
                    Object tagValue = alignment.getAttribute(tag);
                    if (tagValue != null) {
                        c = tagValueColors.get(tagValue.toString());
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

    public static PEStats getPEStats(Alignment alignment, RenderOptions renderOptions) {
        String lb = alignment.getLibrary();
        if (lb == null) lb = "null";
        PEStats peStats = null;
        if (renderOptions.peStats != null) {
            peStats = renderOptions.peStats.get(lb);
        }
        return peStats;
    }

    /**
     * Determine if alignment insert size is outside max or min
     * range for outliers.
     *
     * @param alignment
     * @param renderOptions
     * @return -1 if unknown (stats not computed), 0 if not
     *         an outlier, 1 if outlier
     */
    private int getOutlierStatus(Alignment alignment, RenderOptions renderOptions) {
        PEStats peStats = getPEStats(alignment, renderOptions);
        if (renderOptions.isComputeIsizes() && peStats != null) {
            int minThreshold = peStats.getMinOutlierInsertSize();
            int maxThreshold = peStats.getMaxOutlierInsertSize();
            int dist = Math.abs(alignment.getInferredInsertSize());
            if (dist >= minThreshold || dist <= maxThreshold) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return -1;
        }
    }

    /**
     * Returns -1 if alignment distance is less than minimum,
     * 0 if within bounds, and +1 if above maximum.
     *
     * @param alignment
     * @return
     */
    private int compareToBounds(Alignment alignment, RenderOptions renderOptions) {
        int minThreshold = renderOptions.getMinInsertSize();
        int maxThreshold = renderOptions.getMaxInsertSize();
        PEStats peStats = getPEStats(alignment, renderOptions);
        if (renderOptions.isComputeIsizes() && peStats != null) {
            minThreshold = peStats.getMinThreshold();
            maxThreshold = peStats.getMaxThreshold();
        }

        int dist = Math.abs(alignment.getInferredInsertSize());

        if (dist < minThreshold) return -1;
        if (dist > maxThreshold) return +1;
        return 0;
    }

    /**
     * Assuming we want to color a pair of alignments based on their distance,
     * this returns an appropriate color
     *
     * @param pair
     * @return
     */
    private static Color getColorRelDistance(PairedAlignment pair) {
        if (pair.secondAlignment == null) {
            return grey1;
        }

        int dist = Math.abs(pair.getInferredInsertSize());
        double logDist = Math.log(dist);
        Color minColor = smallISizeColor;
        Color maxColor = largeISizeColor;
        ContinuousColorScale colorScale = new ContinuousColorScale(0, 20, minColor, maxColor);
        return colorScale.getColor((float) logDist);
    }


    /**
     * Returns a color to flag unexpected pair orientations.  Expected orientations (e.g. FR for Illumina) get the
     * neutral grey color
     *
     * @param alignment
     * @param peStats
     * @return
     */
    private Color getOrientationColor(Alignment alignment, PEStats peStats) {
        AlignmentTrack.OrientationType type = getOrientationType(alignment, peStats);
        Color c = typeToColorMap.get(type);
        return c == null ? grey1 : c;

    }

    static AlignmentTrack.OrientationType getOrientationType(Alignment alignment, PEStats peStats) {
        AlignmentTrack.OrientationType type = null;
        if (alignment.isPaired()) { // && !alignment.isProperPair()) {
            final String pairOrientation = alignment.getPairOrientation();
            if (peStats != null) {
                PEStats.Orientation libraryOrientation = peStats.getOrientation();
                switch (libraryOrientation) {
                    case FR:
                        //if (!alignment.isSmallInsert()) {
                        // if the isize < read length the reads might overlap, invalidating this test
                        type = frOrientationTypes.get(pairOrientation);
                        //}
                        break;
                    case RF:
                        type = rfOrientationTypes.get(pairOrientation);
                        break;
                    case F1F2:
                        type = f1f2OrientationTypes.get(pairOrientation);
                        break;
                    case F2F1:
                        type = f2f1OrientationTypes.get(pairOrientation);
                        break;
                }

            } else {
                // No peStats for this library, just guess
                if (alignment.getAttribute("CS") != null) {
                    type = f2f1OrientationTypes.get(pairOrientation);
                } else {
                    type = frOrientationTypes.get(pairOrientation);
                }
            }
        }

        return type;
    }

    /**
     * Similar to "pair orientation" color, but this method does not attempt to interpret orientations.
     *
     * @param alignment
     * @return
     */
    private Color getTemplateStrandColor(Alignment alignment) {

        Color c = null;
        if (alignment.isPaired()) {

            final String pairOrientation = alignment.getPairOrientation();
            return tagValueColors.get(pairOrientation);
        }

        return c == null ? grey1 : c;

    }

    public SortedSet<Shape> curveOverlap(double x) {
        QuadCurve2D tcurve = new QuadCurve2D.Double();
        tcurve.setCurve(x, 0, x, 0, x, 0);
        SortedSet overlap = new TreeSet(arcsByStart.headSet(tcurve, true));
        overlap.retainAll(arcsByEnd.tailSet(tcurve, true));
        return overlap;
    }


    public Alignment getAlignmentForCurve(Shape curve) {
        return curveMap.get(curve);
    }

    public void clearCurveMaps() {
        curveMap.clear();
        arcsByStart.clear();
        arcsByEnd.clear();
    }
}
