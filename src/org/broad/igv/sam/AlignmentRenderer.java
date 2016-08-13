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

package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.MonocolorScale;
import org.broad.igv.sam.AlignmentTrack.ColorOption;
import org.broad.igv.sam.AlignmentTrack.RenderOptions;
import org.broad.igv.sam.AlignmentTrack.ShadeBasesOption;
import org.broad.igv.sam.BisulfiteBaseInfo.DisplayStatus;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.color.*;
import org.broad.igv.util.ChromosomeColors;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.QuadCurve2D;
import java.awt.geom.Rectangle2D;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AlignmentRenderer implements FeatureRenderer {


    public static final HSLColorTable tenXColorTable1 = new HSLColorTable(30);
    public static final HSLColorTable tenXColorTable2 = new HSLColorTable(270);
    public static final GreyscaleColorTable tenXColorTable3 = new GreyscaleColorTable();

    public static final MonocolorScale RED_SCALE = new MonocolorScale(0, 100000, Color.RED);
    public static final MonocolorScale BLUE_SCALE = new MonocolorScale(0, 100000, Color.blue);
    public static final MonocolorScale GRAY_SCALE = new MonocolorScale(0, 100000, Color.GRAY);

    private static Logger log = Logger.getLogger(AlignmentRenderer.class);

    public static final Color GROUP_DIVIDER_COLOR = new Color(200, 200, 200);

    // A "dummy" reference for soft-clipped reads.
    private static byte[] softClippedReference = new byte[1000];

    private static Color smallISizeColor = new Color(0, 0, 150);
    private static Color largeISizeColor = new Color(150, 0, 0);
    private static Color purple = new Color(118, 24, 220);
    private static Color deletionColor = Color.black;
    private static Color skippedColor = new Color(150, 184, 200);
    private static Color unknownGapColor = new Color(0, 150, 0);
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
    private Map<String, ColorTable> tagValueColors;
    private ColorTable defaultTagColors;

    private final Color LR_COLOR = grey1; // "Normal" alignment color
    //private final Color LR_COLOR_12 = new Color(190, 190, 210);
    //private final Color LR_COLOR_21 = new Color(210, 190, 190);
    private final Color RL_COLOR = new Color(0, 150, 0);
    private final Color RR_COLOR = new Color(20, 50, 200);
    private final Color LL_COLOR = new Color(0, 150, 150);

    private final Color OUTLINE_COLOR = new Color(185, 185, 185);
    private static final Color SUPPLEMENTARY_OUTLINE_COLOR = Color.blue;

    private static Map<String, AlignmentTrack.OrientationType> frOrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> f1f2OrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> f2f1OrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> rfOrientationTypes;
    private Map<AlignmentTrack.OrientationType, Color> typeToColorMap;

    public static HashMap<Character, Color> nucleotideColors;

    static {
        initializeTagTypes();
        setNucleotideColors();
    }

    PreferenceManager prefs;

    private static AlignmentRenderer instance;

    private TreeSet<Shape> arcsByStart;
    private TreeSet<Shape> arcsByEnd;
    private HashMap<Shape, Alignment> curveMap;

    private static void setNucleotideColors() {

        PreferenceManager prefs = PreferenceManager.getInstance();

        nucleotideColors = new HashMap();

        Color a = ColorUtilities.stringToColor(prefs.get(PreferenceManager.SAM_COLOR_A), Color.green);
        Color c = ColorUtilities.stringToColor(prefs.get(PreferenceManager.SAM_COLOR_C), Color.blue);
        Color t = ColorUtilities.stringToColor(prefs.get(PreferenceManager.SAM_COLOR_T), Color.red);
        Color g = ColorUtilities.stringToColor(prefs.get(PreferenceManager.SAM_COLOR_G), Color.gray);
        Color n = ColorUtilities.stringToColor(prefs.get(PreferenceManager.SAM_COLOR_N), Color.gray);

        nucleotideColors.put('A', a);
        nucleotideColors.put('a', a);
        nucleotideColors.put('C', c);
        nucleotideColors.put('c', c);
        nucleotideColors.put('T', t);
        nucleotideColors.put('t', t);
        nucleotideColors.put('G', g);
        nucleotideColors.put('g', g);
        nucleotideColors.put('N', n);
        nucleotideColors.put('n', n);

    }

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
        defaultTagColors = new PaletteColorTable(palette);
        tagValueColors = new HashMap();

        typeToColorMap = new HashMap<>(5);
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
                                 Rectangle trackRect,
                                 RenderOptions renderOptions,
                                 boolean leaveMargin,
                                 Map<String, Color> selectedReadNames,
                                 AlignmentCounts alignmentCounts) {

        double origin = context.getOrigin();
        double locScale = context.getScale();
        Font font = FontManager.getFont(10);
        boolean completeReadsOnly = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_COMPLETE_READS_ONLY);

        if ((alignments != null) && (alignments.size() > 0)) {

            int lastPixelDrawn = -1;

            for (Alignment alignment : alignments) {
                // Compute the start and dend of the alignment in pixels
                double pixelStart = ((alignment.getStart() - origin) / locScale);
                double pixelEnd = ((alignment.getEnd() - origin) / locScale);

                // If any any part of the feature fits in the track rectangle draw  it
                if (pixelEnd < rowRect.x || pixelStart > rowRect.getMaxX()) {
                    continue;
                }

                // Optionally only draw alignments that are completely in view
                if (completeReadsOnly) {
                    if (pixelStart < rowRect.x || pixelEnd > rowRect.getMaxX()) {
                        continue;
                    }
                }

                // If the alignment is 3 pixels or less,  draw alignment as a single block,
                // further detail would not be seen and just add to drawing overhead
                // Does the change for Bisulfite kill some machines?
                double pixelWidth = pixelEnd - pixelStart;
                if ((pixelWidth < 4) && !(AlignmentTrack.isBisulfiteColorType(renderOptions.getColorOption()) && (pixelWidth >= 1))) {

                    Color alignmentColor = getAlignmentColor(alignment, renderOptions);

                    // Optimization for really zoomed out views.  If this alignment occupies screen space already taken,
                    // and it is the default color, skip drawing.
                    if (pixelEnd <= lastPixelDrawn && alignmentColor == grey1) {
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
                    drawPairedAlignment((PairedAlignment) alignment, rowRect, trackRect, context, renderOptions, leaveMargin, selectedReadNames, font, alignmentCounts);
                } else if (alignment instanceof LinkedAlignment) {
                    drawLinkedAlignment((LinkedAlignment) alignment, rowRect, trackRect, context, renderOptions, leaveMargin, selectedReadNames, font, alignmentCounts);
                } else {
                    Color alignmentColor = getAlignmentColor(alignment, renderOptions);
                    Graphics2D g = context.getGraphic2DForColor(alignmentColor);
                    g.setFont(font);
                    drawAlignment(alignment, rowRect, trackRect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames, alignmentCounts);
                }
            }

            // Optionally draw a border around the center base
            boolean showCenterLine = prefs.getAsBoolean(PreferenceManager.SAM_SHOW_CENTER_LINE);
            final int bottom = rowRect.y + rowRect.height;
            if (showCenterLine) {
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

    private void drawLinkedAlignment(LinkedAlignment alignment, Rectangle rowRect, Rectangle trackRect, RenderContext context, RenderOptions renderOptions, boolean leaveMargin, Map<String, Color> selectedReadNames, Font font, AlignmentCounts alignmentCounts) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        Color alignmentColor = getAlignmentColor(alignment, renderOptions);

        List<Alignment> barcodedAlignments = alignment.alignments;
        if (barcodedAlignments.size() > 0) {
            Alignment firstAlignment = barcodedAlignments.get(0);
            Graphics2D g = context.getGraphic2DForColor(alignmentColor);
            g.setFont(font);
            if (barcodedAlignments.size() > 1) {
                Color lineColor = new Color(alignmentColor.getRed() / 255f, alignmentColor.getGreen() / 255f, alignmentColor.getBlue() / 255f, 0.3f);
                Graphics2D gline = context.getGraphic2DForColor(lineColor);
                int startX = (int) ((firstAlignment.getEnd() - origin) / locScale);
                int endX = (int) ((barcodedAlignments.get(barcodedAlignments.size() - 1).getStart() - origin) / locScale);
                int h = (int) Math.max(1, rowRect.getHeight() - (leaveMargin ? 2 : 0));
                int y = (int) (rowRect.getY());
                startX = Math.max(rowRect.x, startX);
                endX = Math.min(rowRect.x + rowRect.width, endX);
                gline.drawLine(startX, y + h / 2, endX, y + h / 2);
            }
            for (Alignment al : barcodedAlignments) {
                drawAlignment(al, rowRect, trackRect, g, context, alignmentColor, renderOptions, leaveMargin, selectedReadNames, alignmentCounts);
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

        int d = Math.max(0, (int) (arrowLength + 2 - AlignmentPacker.MIN_ALIGNMENT_SPACING / context.getScale()));
        if (alignment.isNegativeStrand()) x += d;
        if (!alignment.isNegativeStrand()) w -= d;


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
     * @param alignmentCounts
     */
    private void drawPairedAlignment(
            PairedAlignment pair,
            Rectangle rowRect,
            Rectangle trackRect,
            RenderContext context,
            RenderOptions renderOptions,
            boolean leaveMargin,
            Map<String, Color> selectedReadNames,
            Font font,
            AlignmentCounts alignmentCounts) {

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
        drawAlignment(pair.firstAlignment, rowRect, trackRect, g, context, alignmentColor1, renderOptions, leaveMargin, selectedReadNames, alignmentCounts);

        //If the paired alignment is in memory, we draw it.
        //However, we get the coordinates from the first alignment
        if (pair.secondAlignment != null) {

            if (alignmentColor2 == null) {
                alignmentColor2 = getAlignmentColor(pair.secondAlignment, renderOptions);
            }
            g = context.getGraphic2DForColor(alignmentColor2);

            drawAlignment(pair.secondAlignment, rowRect, trackRect, g, context, alignmentColor2, renderOptions, leaveMargin, selectedReadNames, alignmentCounts);
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
     * Draw a single ungapped block in an alignment.
     */
    private void drawAlignmentBlock(Graphics2D blockGraphics, Graphics2D outlineGraphics, Graphics2D terminalGraphics,
                                    boolean isNegativeStrand, int alignmentChromStart, int alignmentChromEnd, int blockChromStart, int blockChromEnd,
                                    int blockPxStart, int blockPxWidth, int y, int h, double locSale) {

        if (blockPxWidth == 0) {
            return;
        } // skip blocks too small to render

        int blockPxEnd = blockPxStart + blockPxWidth;

        boolean leftmost = (blockChromStart == alignmentChromStart),
                rightmost = (blockChromEnd == alignmentChromEnd),
                tallEnoughForArrow = h > 8;


        if (h == 1) {
            blockGraphics.drawLine(blockPxStart, y, blockPxEnd, y);
        } else {
            Shape blockShape;
            int arrowPxWidth = Math.min(5, blockPxWidth / 6);
            int delta = Math.max(0, (int) (arrowPxWidth + 2 - AlignmentPacker.MIN_ALIGNMENT_SPACING / locSale));
            if (leftmost && isNegativeStrand && tallEnoughForArrow) blockPxStart += delta;
            if (rightmost && !isNegativeStrand && tallEnoughForArrow) blockPxEnd -= delta;

            // Draw block as a rectangle; use a pointed hexagon in terminal block to indicate strand.
            int[] xPoly = {blockPxStart - (leftmost && isNegativeStrand && tallEnoughForArrow ? arrowPxWidth : 0),
                    blockPxStart,
                    blockPxEnd,
                    blockPxEnd + (rightmost && !isNegativeStrand && tallEnoughForArrow ? arrowPxWidth : 0),
                    blockPxEnd,
                    blockPxStart},
                    yPoly = {y + h / 2,
                            y,
                            y,
                            y + h / 2,
                            y + h,
                            y + h};
            blockShape = new Polygon(xPoly, yPoly, xPoly.length);

            blockGraphics.fill(blockShape);
            if (outlineGraphics != null) {
                outlineGraphics.draw(blockShape);
            }
        }

        // If the block is too small for a pointed hexagon arrow, then indicate strand with a line.
        if (!tallEnoughForArrow) {
            int tH = Math.max(1, h - 1);
            if (leftmost && isNegativeStrand) {
                terminalGraphics.drawLine(blockPxStart, y, blockPxStart, y + tH);
            }
            if (rightmost && !isNegativeStrand) {
                terminalGraphics.drawLine(blockPxEnd, y, blockPxEnd, y + tH);
            }
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
     * @param alignmentCounts
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
            Map<String, Color> selectedReadNames,
            AlignmentCounts alignmentCounts) {

        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();

        // No blocks.  Note: SAM/BAM alignments always have at least 1 block
        if (blocks == null || blocks.length == 0) {
            drawSimpleAlignment(alignment, rowRect, g, context, renderOptions.flagUnmappedPairs);
            return;
        }

        boolean hideSmallIndelsBP = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_HIDE_SMALL_INDEL_BP);
        int indelThresholdBP = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_SMALL_INDEL_BP_THRESHOLD);
        boolean hideSmallIndelsPixel = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_HIDE_SMALL_INDEL_PIXEL);
        int indexThresholdPixel = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_SMALL_INDELS_PIXEL_THRESHOLD);


        // Scale and position of the alignment rendering.
        double locScale = context.getScale();
        int h = (int) Math.max(1, rowRect.getHeight() - (leaveMargin ? 2 : 0));
        int y = (int) (rowRect.getY());

        // Get a graphics context for outlining alignment blocks.
        Graphics2D outlineGraphics = null;
        if (selectedReadNames.containsKey(alignment.getReadName())) {
            Color c = selectedReadNames.get(alignment.getReadName());
            c = (c == null) ? Color.blue : c;
            outlineGraphics = context.getGraphic2DForColor(c);
            outlineGraphics.setStroke(thickStroke);
        } else if (renderOptions.flagUnmappedPairs && alignment.isPaired() && !alignment.getMate().isMapped()) {
            outlineGraphics = context.getGraphic2DForColor(Color.red);
        } else if (alignment.getMappingQuality() == 0 && renderOptions.flagZeroQualityAlignments) {
            outlineGraphics = context.getGraphic2DForColor(OUTLINE_COLOR);
        }

        // Define a graphics context for small insertions.
        Graphics2D smallInsertionGraphics = context.getGraphic2DForColor(purple);

        // Define a graphics context for indel labels.
        Graphics2D largeIndelGraphics = (Graphics2D) context.getGraphics().create();
        largeIndelGraphics.setFont(FontManager.getFont(Font.BOLD, h - 2));
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
            largeIndelGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }

        // Get a graphics context for drawing individual basepairs.
        Graphics2D bpGraphics = context.getGraphic2D("BASE");
        int dX = (int) Math.max(1, (1.0 / locScale));
        if (dX >= 8) {
            Font f = FontManager.getFont(Font.BOLD, Math.min(dX, 12));
            bpGraphics.setFont(f);
        }

        // Define the graphics contexts for various types of gap.
        Graphics2D defaultGapGraphics = context.getGraphic2DForColor(deletionColor);
        defaultGapGraphics.setStroke(thickStroke);
        Graphics2D unknownGapGraphics = context.getGraphic2DForColor(unknownGapColor);
        Graphics2D skippedRegionGapGraphics = context.getGraphic2DForColor(skippedColor);

        // Get a graphics context to indicate the end of a read.
        Graphics2D terminalGraphics = context.getGraphic2DForColor(Color.DARK_GRAY);

        /* Process the alignment. */
        AlignmentBlock firstBlock = blocks[0], lastBlock = blocks[blocks.length - 1];
        int alignmentChromStart = (int) firstBlock.getStart(),
                alignmentChromEnd = (int) (lastBlock.getStart() + lastBlock.getLength());

        // BED-style coordinate for the visible context.  Do not draw outside the context.
        int contextChromStart = (int) context.getOrigin(),
                contextChromEnd = (int) context.getEndLocation();
        // BED-style start coordinate for the next alignment block to draw.
        int blockChromStart = (int) Math.max(alignmentChromStart, contextChromStart);

        // Draw aligment blocks separated by gaps.  Define the blocks by walking through the gap list,
        // skipping over gaps that are too small to show at the curren resolution.
        java.util.List<Gap> gaps = alignment.getGaps();
        if (gaps != null) {

            for (Gap gap : gaps) {
                int gapChromStart = (int) gap.getStart(),
                        gapChromWidth = (int) gap.getnBases(),
                        gapChromEnd = gapChromStart + gapChromWidth,
                        gapPxWidth = (int) Math.max(1, gapChromWidth / locScale),
                        gapPxEnd = (int) ((Math.min(contextChromEnd, gapChromEnd) - contextChromStart) / locScale);

                if (gapChromEnd <= contextChromStart) { // gap ends before the visible context
                    continue; // move to next gap
                } else if (gapChromStart >= contextChromEnd) { // gap starts after the visible context
                    break; // done examining gaps
                }

                // Draw the gap if it is sufficiently large at the current zoom.
                boolean drawGap = ((!hideSmallIndelsBP || gapChromWidth > indelThresholdBP) &&
                        (!hideSmallIndelsPixel || gapPxWidth >= indexThresholdPixel));
                if (!drawGap) {
                    continue;
                }

                // Draw the preceding alignment block.
                int blockPxStart = (int) ((blockChromStart - contextChromStart) / locScale),
                        blockChromEnd = gapChromStart,
                        blockPxWidth = (int) Math.max(1, (blockChromEnd - blockChromStart) / locScale),
                        blockPxEnd = blockPxStart + blockPxWidth;
                drawAlignmentBlock(g, outlineGraphics, terminalGraphics, alignment.isNegativeStrand(),
                        alignmentChromStart, alignmentChromEnd, blockChromStart, blockChromEnd,
                        blockPxStart, blockPxWidth, y, h, locScale);

                // Draw the gap line.
                Graphics2D gapGraphics = defaultGapGraphics;
                if (gap.getType() == SAMAlignment.UNKNOWN) {
                    gapGraphics = unknownGapGraphics;
                } else if (gap.getType() == SAMAlignment.SKIPPED_REGION) {
                    gapGraphics = skippedRegionGapGraphics;
                }

                gapGraphics.drawLine(blockPxEnd + 1, y + h / 2, gapPxEnd, y + h / 2);

                // Label the size of the deletion if it is "large" and the label fits.
                if (renderOptions.isFlagLargeIndels() && gapChromWidth > renderOptions.getLargeInsertionsThreshold()) {
                    drawLargeIndelLabel(largeIndelGraphics, false, Globals.DECIMAL_FORMAT.format(gapChromWidth), (int) ((blockPxEnd + gapPxEnd) / 2), y, h, gapPxEnd - blockPxEnd - 2);
                }

                // Start the next alignment block after the gap.
                blockChromStart = gapChromEnd;
            }
        }

        // Draw the final block after the last gap.
        int blockPxStart = (int) ((blockChromStart - contextChromStart) / locScale),
                blockChromEnd = (int) Math.min(contextChromEnd, alignmentChromEnd),
                blockPxWidth = (int) Math.max(1, (blockChromEnd - blockChromStart) / locScale),
                blockPxEnd = blockPxStart + blockPxWidth,
                lastBlockPxEnd = blockPxEnd;
        drawAlignmentBlock(g, outlineGraphics, terminalGraphics, alignment.isNegativeStrand(),
                alignmentChromStart, alignmentChromEnd, blockChromStart, blockChromEnd,
                blockPxStart, blockPxWidth, y, h, locScale);

        // Draw insertions.
        drawInsertions(contextChromStart, rowRect, locScale, alignment, smallInsertionGraphics, largeIndelGraphics, renderOptions);

        // Draw basepairs / mismatches.
        if (locScale < 100) {
            if (renderOptions.showMismatches || renderOptions.showAllBases) {
                boolean quickConsensus = prefs.getAsBoolean(PreferenceManager.SAM_QUICK_CONSENSUS_MODE);
                for (AlignmentBlock aBlock : alignment.getAlignmentBlocks()) {
                    int aBlockChromStart = (int) aBlock.getStart(),
                            aBlockChromEnd = (int) (aBlock.getStart() + aBlock.getLength());

                    if (aBlockChromEnd <= contextChromStart) { // block ends before the visible context
                        continue; // move to next block
                    } else if (aBlockChromStart >= contextChromEnd) { // block starts after the visible context
                        break; // done examining blocks
                    }

                    drawBases(context, bpGraphics, rowRect, alignment, aBlock, alignmentCounts, quickConsensus, alignmentColor, renderOptions);
                }
            }
        }

        //Draw straight line up for viewing arc pairs, if mate on a different chromosome
        if (renderOptions.isPairedArcView()) {
            try {
                Graphics2D gLine = context.getGraphic2DForColor(alignmentColor);
                if (!alignment.getChr().equalsIgnoreCase(alignment.getMate().getChr())) {
                    gLine.drawLine(lastBlockPxEnd, y + h / 2, lastBlockPxEnd, (int) trackRect.getMinY());
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
     * @param g
     * @param rect
     * @param baseAlignment
     * @param block
     * @param alignmentCounts
     * @param quickConsensus
     * @param alignmentColor
     * @param renderOptions
     */
    private void drawBases(RenderContext context,
                           Graphics2D g,
                           Rectangle rect,
                           Alignment baseAlignment,
                           AlignmentBlock block,
                           AlignmentCounts alignmentCounts,
                           boolean quickConsensus,
                           Color alignmentColor,
                           RenderOptions renderOptions) {

        boolean isSoftClipped = block.isSoftClipped();


        // Get the base qualities, start/end,  and reference sequence
        String chr = context.getChr();
        final int start = block.getStart();
        final int end = block.getEnd();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        final byte[] reference = isSoftClipped ? softClippedReference : genome.getSequence(chr, start, end);

        boolean haveBases = (block.hasBases() && block.getLength() > 0);

        ShadeBasesOption shadeBasesOption = renderOptions.shadeBasesOption;
        ColorOption colorOption = renderOptions.getColorOption();

        // Disable showAllBases in bisulfite mode
        boolean showAllBases = renderOptions.showAllBases &&
                !(colorOption == ColorOption.BISULFITE || colorOption == ColorOption.NOMESEQ);

        if (!showAllBases && (!haveBases || reference == null)) {
            return;
        }

        byte[] read;
        if (haveBases) {
            read = block.getBases();
        } else {
            read = reference;
        }


        double locScale = context.getScale();
        double origin = context.getOrigin();

        // Compute bounds
        int pY = (int) rect.getY();
        int dY = (int) rect.getHeight();
        int dX = (int) Math.max(1, (1.0 / locScale));

        BisulfiteBaseInfo bisinfo = null;
        boolean nomeseqMode = (renderOptions.getColorOption().equals(AlignmentTrack.ColorOption.NOMESEQ));
        boolean bisulfiteMode = AlignmentTrack.isBisulfiteColorType(renderOptions.getColorOption());
        if (nomeseqMode) {
            bisinfo = new BisulfiteBaseInfoNOMeseq(reference, baseAlignment, block, renderOptions.bisulfiteContext);
        } else if (bisulfiteMode) {
            bisinfo = new BisulfiteBaseInfo(reference, baseAlignment, block, renderOptions.bisulfiteContext);
        }


        for (int loc = start; loc < end; loc++) {
            int idx = loc - start;

            boolean misMatch = haveBases && AlignmentUtils.isMisMatch(reference, read, isSoftClipped, idx);

            if (showAllBases || (!bisulfiteMode && misMatch) ||
                    (bisulfiteMode && (!DisplayStatus.NOTHING.equals(bisinfo.getDisplayStatus(idx))))) {
                char c = (char) read[idx];

                Color color = nucleotideColors.get(c);
                if (bisulfiteMode) color = bisinfo.getDisplayColor(idx);
                if (color == null) {
                    color = Color.black;
                }

                if (ShadeBasesOption.QUALITY == shadeBasesOption) {
                    byte qual = block.getQuality(loc - start);
                    color = getShadedColor(qual, color, alignmentColor, prefs);
                } else if (ShadeBasesOption.FLOW_SIGNAL_DEVIATION_READ == shadeBasesOption || ShadeBasesOption.FLOW_SIGNAL_DEVIATION_REFERENCE == shadeBasesOption) {
                    if (block.hasFlowSignals()) {
                        color = getFlowSignalColor(reference, misMatch, genome, block, chr, start, loc, idx, shadeBasesOption, alignmentColor, color);
                    }
                }

                double bisulfiteXaxisShift = (bisulfiteMode) ? bisinfo.getXaxisShift(idx) : 0;

                // If there is room for text draw the character, otherwise
                // just draw a rectangle to represent the
                int pX = (int) (((double) loc + bisulfiteXaxisShift - origin) / locScale);

                // Don't draw out of clipping rect
                if (pX > rect.getMaxX()) {
                    break;
                } else if (pX + dX < rect.getX()) {
                    continue;
                }

                BisulfiteBaseInfo.DisplayStatus bisstatus = (bisinfo == null) ? null : bisinfo.getDisplayStatus(idx);
                if (isSoftClipped || bisulfiteMode ||
                        // In "quick consensus" mode, only show mismatches at positions with a consistent alternative basepair.
                        (!quickConsensus || alignmentCounts.isMismatch(loc, reference[idx], chr, prefs.getAsFloat(PreferenceManager.SAM_ALLELE_THRESHOLD)))
                        ) {
                    drawBase(g, color, c, pX, pY, dX, dY, bisulfiteMode, bisstatus);
                }
            }
        }

    }

    /**
     * NB: this may estimate the reference homopolymer length incorrect in some cases, especially when we have
     * an overcall/undercall situation.  Proper estimation of the reads observed versus expected homopolymer
     * length should use flow signal alignment (SamToFlowspace): https://github.com/iontorrent/Ion-Variant-Hunter
     *
     * @param reference
     * @param misMatch
     * @param genome
     * @param block
     * @param chr
     * @param start
     * @param loc
     * @param idx
     * @param shadeBasesOption
     * @param alignmentColor
     * @param color
     * @return
     */
    private Color getFlowSignalColor(byte[] reference, boolean misMatch, Genome genome,
                                     AlignmentBlock block, String chr, int start, int loc, int idx,
                                     ShadeBasesOption shadeBasesOption,
                                     Color alignmentColor, Color color) {
        int flowSignal = (int) block.getFlowSignalSubContext(loc - start).getCurrentSignal();
        int expectedFlowSignal;
        if (ShadeBasesOption.FLOW_SIGNAL_DEVIATION_READ == shadeBasesOption) {
            expectedFlowSignal = 100 * (short) ((flowSignal + 50.0) / 100.0);
        } else {
            // NB: this may estimate the reference homopolymer length incorrect in some cases, especially when we have
            // an overcall/undercall situation.  Proper estimation of the reads observed versus expected homopolymer
            // length should use flow signal alignment (SamToFlowspace): https://github.com/iontorrent/Ion-Variant-Hunter
            if (!misMatch) {
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
        if (Byte.MAX_VALUE < flowSignalDiff) {
            qual = Byte.MAX_VALUE;
        } else if (flowSignalDiff < Byte.MIN_VALUE) {
            qual = Byte.MIN_VALUE;
        } else {
            qual = (byte) flowSignalDiff;
        }
        // Finally, get the color
        return getShadedColor(qual, color, alignmentColor, prefs);
    }

    /**
     * Draw the base using either a letter or character, using the given color,
     * depending on size and bisulfite status
     *
     * @param g
     * @param color
     * @param c
     * @param pX
     * @param pY
     * @param dX
     * @param dY
     * @param bisulfiteMode
     * @param bisstatus
     */
    private void drawBase(Graphics2D g, Color color, char c, int pX, int pY, int dX, int dY, boolean bisulfiteMode,
                          DisplayStatus bisstatus) {
        if (((dY >= 12) && (dX >= 8)) && (!bisulfiteMode || (bisulfiteMode && bisstatus.equals(DisplayStatus.CHARACTER)))) {
            g.setColor(color);
            GraphicUtils.drawCenteredText(new char[]{c}, pX, pY + 1, dX, dY - 2, g);
        } else {

            int pX0i = pX, dXi = dX;

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

    private void drawLargeIndelLabel(Graphics2D g, boolean isInsertion, String labelText, int pxCenter, int pxTop, int pxH, int pxWmax) {
        final int pxPad = 2;   // text padding in the label
        final int pxWing = 2;  // width of the cursor "wing"

        // Calculate the width required to draw the label
        Rectangle2D textBounds = g.getFontMetrics().getStringBounds(labelText, g);
        int pxTextW = 2 * pxPad + (int) textBounds.getWidth();
        boolean doesTextFit = (pxTextW < pxWmax);

        if (!doesTextFit && !isInsertion) {
            return;
        } // only label deletions when the text fits

        // Calculate the pixel bounds of the label
        int pxW = (int) Math.max(2, Math.min(pxTextW, pxWmax)),
                pxLeft = pxCenter - (int) Math.ceil(pxW / 2),
                pxRight = pxLeft + pxW;

        // Draw the label
        g.setColor(isInsertion ? purple : Color.white);
        g.fillRect(pxLeft, pxTop, pxRight - pxLeft, pxH);

        // TODO -- record this "object" for popup text
        if (isInsertion) {
            g.fillRect(pxLeft - pxWing, pxTop, pxRight - pxLeft + 2 * pxWing, 2);
            g.fillRect(pxLeft - pxWing, pxTop + pxH - 2, pxRight - pxLeft + 2 * pxWing, 2);
        } // draw "wings" For insertions

        if (doesTextFit) {
            g.setColor(isInsertion ? Color.white : purple);
            g.drawString(labelText, pxLeft + pxPad, pxTop + pxH - 2);
        } // draw the text if it fits
    }

    private void drawInsertions(double origin, Rectangle rect, double locScale, Alignment alignment,
                                Graphics2D gSmallInsertion, Graphics2D gLargeInsertion, RenderOptions renderOptions) {

        AlignmentBlock[] insertions = alignment.getInsertions();
        if (insertions != null) {

            boolean hideSmallIndelsBP = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_HIDE_SMALL_INDEL_BP);
            int indelThresholdBP = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_SMALL_INDEL_BP_THRESHOLD);
            boolean hideSmallIndelsPixel = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_HIDE_SMALL_INDEL_PIXEL);
            int indexThresholdPixel = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_SMALL_INDELS_PIXEL_THRESHOLD);


            for (AlignmentBlock aBlock : insertions) {

                int x = (int) ((aBlock.getStart() - origin) / locScale);
                int pxWidth = (int) (aBlock.getBases().length / locScale);
                int h = (int) Math.max(1, rect.getHeight() - 2);
                int y = (int) (rect.getY() + (rect.getHeight() - h) / 2) - 1;

                // Don't draw out of clipping rect
                if (x > rect.getMaxX()) {
                    break;
                } else if (x < rect.getX()) {
                    continue;
                }

                if ((!hideSmallIndelsBP || aBlock.getBases().length > indelThresholdBP) &&
                        (!hideSmallIndelsPixel || pxWidth >= indexThresholdPixel)) {
                    if (renderOptions.isFlagLargeIndels() && aBlock.getBases().length > renderOptions.getLargeInsertionsThreshold()) {
                        drawLargeIndelLabel(gLargeInsertion, true, Globals.DECIMAL_FORMAT.format(aBlock.getBases().length), x - 1, y, h, pxWidth);
                    } else {
                        gSmallInsertion.fillRect(x - 2, y, 4, 2);
                        gSmallInsertion.fillRect(x - 1, y, 2, h);
                        gSmallInsertion.fillRect(x - 2, y + h - 2, 4, 2);
                    }
                }
            }
        }
    }

    private Color getAlignmentColor(Alignment alignment, RenderOptions renderOptions) {

        // Set color used to draw the feature.  Highlight features that intersect the
        // center line.  Also restorePersistentState row "score" if alignment intersects center line


        Color color = alignment.getColor();
        if (color != null) return color;   // Color has been explicitly set

        Color c = grey1;

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
                if (colorOption == ColorOption.PAIR_ORIENTATION) {
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

                        ColorTable ctable;
                        String ctableKey;

                        String groupByTag = renderOptions.getGroupByTag();
                        if (groupByTag == null) {
                            ctable = defaultTagColors;

                        } else {
                            Object g = alignment.getAttribute(groupByTag);
                            String group = g == null ? "" : g.toString();
                            ctableKey = groupByTag + ":" + group;
                            ctable = tagValueColors.get(ctableKey);
                            if (ctable == null) {

                                if (groupByTag.equals("HP")) {
                                    ctable = getTenXColorTable(group);
                                } else {
                                    ctable = defaultTagColors;
                                }

                                tagValueColors.put(group, ctable);
                            }
                        }

                        c = ctable.get(tagValue.toString());

                    }
                }
                break;
            case LINK_STRAND:
                if (alignment instanceof LinkedAlignment && ((LinkedAlignment) alignment).getStrand() == Strand.NONE) {
                    c = LL_COLOR;
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

    private ColorTable getTenXColorTable(String group) {
        ColorTable ctable;
        if (group.equals("1")) {
            ctable = tenXColorTable1;

        } else if (group.equals("2")) {
            ctable = tenXColorTable2;
        } else {
            ctable = tenXColorTable3;
        }
        return ctable;
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
     * an outlier, 1 if outlier
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
