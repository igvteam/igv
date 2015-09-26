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

package org.broad.igv.renderer;


import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.collections.MultiMap;
import org.broad.igv.variant.VariantRenderer;

import java.awt.*;
import java.awt.font.LineMetrics;
import java.awt.geom.Rectangle2D;
import java.util.*;
import java.util.List;


/**
 * Renderer for the full "IGV" feature interface
 */
public class IGVFeatureRenderer extends FeatureRenderer {

    private static Logger log = Logger.getLogger(IGVFeatureRenderer.class);

    // Constants
    static protected final int NORMAL_STRAND_Y_OFFSET = 14;
    static protected final int ARROW_SPACING = 30;
    static protected final int NO_STRAND_THICKNESS = 2;
    static protected final int REGION_STRAND_THICKNESS = 4;
    static final int BLOCK_HEIGHT = 14;
    static final int THIN_BLOCK_HEIGHT = 6;
    protected Color AA_COLOR_1 = new Color(92, 92, 164);
    protected Color AA_COLOR_2 = new Color(12, 12, 120);
    public static final Color DULL_BLUE = new Color(0, 0, 200);
    public static final Color DULL_RED = new Color(200, 0, 0);
    static final int NON_CODING_HEIGHT = 8;

    private static final Color VARIANT_HET_COLOR = Color.blue.brighter();
    private static final Color VARIANT_HOM_COLOR = new Color(0, 245, 255);

    float viewLimitMin = Float.NaN;
    float viewLimitMax = Float.NaN;

    // Use the max of these values to determine where
    // text should be drawn
    //protected double lastFeatureLineMaxY = 0;
    //protected double lastFeatureBoundsMaxY = 0;
    //protected double lastRegionMaxY = 0;

    protected boolean drawBoundary = false;

    int blockHeight = BLOCK_HEIGHT;
    int thinBlockHeight = THIN_BLOCK_HEIGHT;

    //Map from Exon to y offset
    //Could use more coordinates, but they are all contained in the Exon
    private Map<IExon, Integer> exonMap = new HashMap<IExon, Integer>(100);
    private Set<String> drawnNames = new HashSet<String>(100);

    protected boolean isGenotypeRenderer = false;

    private static final int MAX_NAME_LENGTH = 60;


    /**
     * Note:  assumption is that featureList is sorted by start position.
     *
     * @param featureList
     * @param context
     * @param trackRectangle
     * @param track
     */
    public void render(List<IGVFeature> featureList,
                       RenderContext context,
                       Rectangle trackRectangle,
                       Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();
        double end = origin + trackRectangle.getWidth() * locScale;

        final Track.DisplayMode displayMode = track.getDisplayMode();
        blockHeight = displayMode == Track.DisplayMode.SQUISHED ? BLOCK_HEIGHT / 2 : BLOCK_HEIGHT;
        thinBlockHeight = displayMode == Track.DisplayMode.SQUISHED ? THIN_BLOCK_HEIGHT / 2 : THIN_BLOCK_HEIGHT;

        // TODO -- use enum instead of string "Color"
        if ((featureList != null) && (featureList.size() > 0)) {

            // Create a graphics object to draw font names.  Graphics are not cached
            // by font, only by color, so its neccessary to create a new one to prevent
            // affecting other tracks.
            Font font = FontManager.getFont(track.getFontSize());
            Graphics2D fontGraphics = (Graphics2D) context.getGraphic2DForColor(Color.BLACK).create();

            if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
                fontGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            }


            fontGraphics.setFont(font);

            // Track coordinates
            double trackRectangleX = trackRectangle.getX();
            double trackRectangleMaxX = trackRectangle.getMaxX();
            double trackRectangleY = trackRectangle.getY();
            double trackRectangleMaxY = trackRectangle.getMaxY();

            int lastNamePixelEnd = -9999;
            int lastPixelEnd = -1;
            int occludedCount = 0;
            int maxOcclusions = 2;


            boolean alternateExonColor = (track instanceof FeatureTrack && ((FeatureTrack) track).isAlternateExonColor());

            for (IGVFeature feature : featureList) {

                if (feature.getEnd() < origin) continue;
                if (feature.getStart() > end) break;

                // Get the pStart and pEnd of the entire feature
                // At extreme zoom levels the
                // virtual pixel value can be too large for an int, so the computation is
                // done in double precision and cast to an int only when its confirmed its
                // within the field of view.
                double virtualPixelStart = (feature.getStart() - origin) / locScale;
                double virtualPixelEnd = (feature.getEnd() - origin) / locScale;

                int pixelStart = (int) Math.round(Math.max(trackRectangleX, virtualPixelStart));
                int pixelEnd = (int) Math.round(Math.min(trackRectangleMaxX, virtualPixelEnd));

                final int pixelWidth = pixelEnd - pixelStart;
                if (isGenotypeRenderer) {
                    if (pixelWidth < 3) {
                        double dx = 3.0 - pixelWidth;
                        pixelStart -= dx / 2.0;
                        pixelEnd += dx / 2.0;
                    }
                }

                // Draw a maximum of "maxOcclusion" small features on top of each other.
                if (pixelEnd <= lastPixelEnd) {
                    if (occludedCount >= maxOcclusions && (pixelEnd - pixelStart) < 3) {
                        continue;
                    } else {
                        occludedCount++;
                    }
                } else {
                    occludedCount = 0;
                    lastPixelEnd = pixelEnd;
                }

                if (isGenotypeRenderer) {
                    renderGenotypeFeature(context, feature, trackRectangle, pixelStart, pixelEnd);
                    continue;
                }

                Color color = getFeatureColor(feature, track);
                Graphics2D g2D = context.getGraphic2DForColor(color);

                // Draw block representing entire feature
                int pixelThickStart = pixelStart;
                int pixelThickEnd = pixelEnd;
                boolean hasExons = false;
                if (feature instanceof BasicFeature) {
                    BasicFeature bf = (BasicFeature) feature;
                    pixelThickStart = (int) Math.max(trackRectangleX, Math.round((bf.getThickStart() - origin) / locScale));
                    pixelThickEnd = (int) Math.min(trackRectangleMaxX, Math.round((bf.getThickEnd() - origin) / locScale));
                    hasExons = bf.hasExons();
                }

                 // Add directional arrows and exons, if there is room.
                int pixelYCenter = trackRectangle.y + NORMAL_STRAND_Y_OFFSET / 2;

                if (hasExons) {
                    if ((pixelWidth < 3)) {
                        drawFeatureBounds(pixelStart, pixelEnd, pixelYCenter, g2D);
                    } else {
                        drawExons(feature, pixelYCenter, context, g2D, trackRectangle, displayMode,
                                alternateExonColor, track.getColor(), track.getAltColor());
                    }
                } else {

                    drawFeatureBlock(pixelStart, pixelEnd, pixelThickStart, pixelThickEnd, pixelYCenter, g2D);
                    Graphics2D arrowGraphics = context.getGraphic2DForColor(Color.WHITE);
                    drawStrandArrows(feature.getStrand(), pixelStart, pixelEnd, pixelYCenter, 0,
                            displayMode, trackRectangle, arrowGraphics);

                    // This is ugly, but alternatives are probably worse
                    if (feature instanceof EncodePeakFeature && pixelWidth > 5) {
                        int peakPosition = ((EncodePeakFeature) feature).getPeakPosition();
                        if (peakPosition > 0) {
                            Color c = g2D.getColor();
                            int peakPixelPosition = (int) ((peakPosition - origin) / locScale);
                            Color peakColor = c == Color.black ? Color.red : Color.black;
                            //if (track.isUseScore()) {
                            //    float alpha = getAlpha(0, 500, feature.getScore());
                            //    peakColor = ColorUtilities.getCompositeColor(DULL_RED, alpha);
                            //}
                            g2D.setColor(peakColor);
                            int pw = Math.min(3, pixelWidth / 5);
                            g2D.fillRect(peakPixelPosition - pw/2, pixelYCenter - thinBlockHeight/2 - 1, pw, thinBlockHeight + 2);
                            g2D.setColor(c);
                        }
                    }

                }


                // Draw name , if there is room
                if (displayMode != Track.DisplayMode.SQUISHED) {
                    String name = feature.getName();
                    if (name != null) {
                        // Limit name display length
                        if (name.length() >= MAX_NAME_LENGTH) {
                            name = name.substring(0, MAX_NAME_LENGTH - 3) + " ...";
                        }

                        LineMetrics lineMetrics = font.getLineMetrics(name, g2D.getFontRenderContext());
                        int fontHeight = (int) Math.ceil(lineMetrics.getHeight());


                        // Draw feature name.  Center it over the feature extent,
                        // or if the feature extends beyond the track bounds over
                        // the track rectangle.
                        int nameStart = Math.max(0, pixelStart);
                        int nameEnd = Math.min(pixelEnd, (int) trackRectangle.getWidth());
                        int textBaselineY = trackRectangle.y + trackRectangle.height - 3;

                        // Calculate the minimum amount of vertical track
                        // space required be we  draw the
                        // track name without drawing over the features
                        int verticalSpaceRequiredForText = textBaselineY - (int) trackRectangleY;

                        if (verticalSpaceRequiredForText <= trackRectangle.height) {
                            lastNamePixelEnd = drawFeatureName(feature, track.getDisplayMode(), nameStart, nameEnd,
                                    lastNamePixelEnd, fontGraphics, textBaselineY);
                        }
                    }
                }

                // If this is the highlight feature highlight it
                if (getHighlightFeature() == feature) {
                    int yStart = pixelYCenter - blockHeight / 2 - 1;
                    Graphics2D highlightGraphics = context.getGraphic2DForColor(Color.cyan);
                    highlightGraphics.drawRect(pixelStart - 1, yStart, (pixelWidth + 2), blockHeight + 2);


                }


            }

            if (drawBoundary) {
                Graphics2D g2D = context.getGraphic2DForColor(Color.LIGHT_GRAY);
                g2D.drawLine((int) trackRectangleX, (int) trackRectangleMaxY - 1,
                        (int) trackRectangleMaxX, (int) trackRectangleMaxY - 1);
            }

            fontGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        }
    }

    protected void renderGenotypeFeature(RenderContext context, IGVFeature feature, Rectangle trackRectangle, int pixelStart, int pixelEnd) {
        Color color = feature.getName().equals("HET") ? VARIANT_HET_COLOR : VARIANT_HOM_COLOR;
        Graphics2D g2D = context.getGraphic2DForColor(color);
        g2D.fillRect(pixelStart, trackRectangle.y, (pixelEnd - pixelStart), trackRectangle.height);
    }

    /**
     * @param pixelStart
     * @param pixelEnd
     * @param pixelThickStart
     * @param pixelThickEnd
     * @param yOffset
     * @param g
     */
    private void drawFeatureBlock(int pixelStart, int pixelEnd, int pixelThickStart, int pixelThickEnd,
                                  int yOffset, Graphics2D g) {

        Graphics2D g2D = (Graphics2D) g.create();

        if (pixelThickStart > pixelStart) {
            g2D.fillRect(pixelStart, yOffset - (thinBlockHeight) / 2,
                    Math.max(1, pixelThickStart - pixelStart), (thinBlockHeight));
        }
        if (pixelThickEnd > 0 && pixelThickEnd < pixelEnd) {
            g2D.fillRect(pixelThickEnd, yOffset - (thinBlockHeight) / 2,
                    Math.max(1, pixelEnd - pixelThickEnd), (thinBlockHeight));
        }

        g2D.fillRect(pixelThickStart, yOffset - (blockHeight - 4) / 2,
                Math.max(1, pixelThickEnd - pixelThickStart), (blockHeight - 4));

        g2D.dispose();
    }

    final private void drawConnectingLine(int startX, int startY, int endX, int endY, Strand strand, Graphics2D g) {

        Graphics2D g2D = (Graphics2D) g.create();
        if (strand == null) {

            BasicStroke befStroke = (BasicStroke) g.getStroke();
            float lineThickness = befStroke.getLineWidth();

            // Double the line thickness
            lineThickness *= NO_STRAND_THICKNESS;
            Stroke stroke = new BasicStroke(lineThickness);
            g2D.setStroke(stroke);
        }
        g2D.drawLine(startX, startY, endX, endY);
        g2D.dispose();

    }


    // If the width is < 3 there isn't room to draw the
    // feature, or orientation.  If the feature has any exons
    // at all indicate by filling a small rect

    private void drawFeatureBounds(int pixelStart, int pixelEnd, int yOffset, Graphics2D g2D) {

        if (pixelEnd == pixelStart) {
            int yStart = yOffset - blockHeight / 2;
            g2D.drawLine(pixelStart, yStart, pixelStart, yStart + blockHeight);
        } else {
            g2D.fillRect(pixelStart, yOffset - blockHeight / 2, pixelEnd - pixelStart, blockHeight);
        }

    }

    protected void drawExons(IGVFeature gene, int yOffset, RenderContext context,
                             Graphics2D g2D, Rectangle trackRectangle, Track.DisplayMode mode,
                             boolean alternateExonColor, Color color1, Color color2) {

        Graphics2D exonNumberGraphics = (Graphics2D) g2D.create();


        exonNumberGraphics.setColor(Color.BLACK);
        exonNumberGraphics.setFont(FontManager.getFont(Font.BOLD, 8));

        // Now get the individual regions of the
        // feature are drawn here

        double theOrigin = context.getOrigin();
        double locationScale = context.getScale();

        Graphics2D fontGraphics = context.getGraphic2DForColor(Color.WHITE);

        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
            exonNumberGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            fontGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }

        boolean colorToggle = true;

        /**
         * We draw connecting lines exon-by-exon,
         * need to keep track of the previous start/end locations
         */
        int lastExonEndX = Integer.MIN_VALUE;
        int lastY = Integer.MIN_VALUE;
        int maxLineEndX = Integer.MIN_VALUE;
        IExon lastExon = null;

        int exonCount = gene.getExons().size();

        for (int idx = 0; idx < exonCount; idx++) {

            Exon exon = gene.getExons().get(idx);

            // Parse expression from tags, if available
            //Credit Michael Poidinger and Solomonraj Wilson, Singapore Immunology Network.
            Float exprValue = null;

            MultiMap<String, String> attributes = exon.getAttributes();
            if (attributes != null && attributes.containsKey("expr")) {
                try {
                    exprValue = Float.parseFloat(attributes.get("expr"));
                } catch (NumberFormatException e) {
                    log.error("Error parsing expression tag " + attributes.get("expr"), e);
                }
            }

            if (exprValue != null) {
                ContinuousColorScale colorScale = PreferenceManager.getInstance().getColorScale(TrackType.GENE_EXPRESSION);
                Color chartColor = colorScale.getColor(exprValue);
                g2D = context.getGraphic2DForColor(chartColor);
            }

            // Added by Solomon - End

            Graphics2D blockGraphics = g2D;
            Graphics2D edgeGraphics = context.getGraphic2DForColor(Color.gray);
            if (alternateExonColor) {
                Color color = colorToggle ? color1 : color2;
                blockGraphics = context.getGraphic2DForColor(color);
                colorToggle = !colorToggle;
            }

            int curYOffset = yOffset;
            boolean drawConnectingLine = true;

            int pStart = getPixelFromChromosomeLocation(exon.getChr(), exon.getStart(), theOrigin, locationScale);
            int pEnd = getPixelFromChromosomeLocation(exon.getChr(), exon.getEnd(), theOrigin, locationScale);


            Graphics2D arrowGraphics = context.getGraphic2DForColor(Color.blue);
            //We draw connecting lines/arrows from previous exon to this one.
            //Exons may not be strictly sorted, we avoid double-drawing using maxLineEndX
            if (drawConnectingLine && lastExonEndX > Integer.MIN_VALUE && lastY > Integer.MIN_VALUE
                    && lastExonEndX >= maxLineEndX) {
                drawConnectingLine(lastExonEndX, lastY, pStart, curYOffset, exon.getStrand(), blockGraphics);
                double angle = Math.atan(-(curYOffset - lastY) / ((pStart - lastExonEndX) + 1e-12));
                drawStrandArrows(gene.getStrand(), lastExonEndX, pStart, lastY, angle, mode, trackRectangle, arrowGraphics);
                maxLineEndX = Math.max(maxLineEndX, pStart);
            }
            lastExonEndX = pEnd;
            lastY = curYOffset;

            if ((pEnd >= trackRectangle.getX()) && (pStart <= trackRectangle.getMaxX())) {
                int pCdStart =
                        Math.min(pEnd, Math.max(pStart,
                                getPixelFromChromosomeLocation(exon.getChr(),
                                        exon.getCdStart(), theOrigin, locationScale)));
                int pCdEnd = Math.max(pStart, Math.min(pEnd,
                        getPixelFromChromosomeLocation(exon.getChr(),
                                exon.getCdEnd(), theOrigin, locationScale)));

                // Entire exon is UTR
                if (exon.isNonCoding()) {
                    int pClippedStart = (int) Math.max(pStart, trackRectangle.getX());
                    int pClippedEnd = (int) Math.min(pEnd, trackRectangle.getMaxX());
                    int pClippedWidth = pClippedEnd - pClippedStart;
                    drawExonRect(blockGraphics, exon, pClippedStart, curYOffset - NON_CODING_HEIGHT / 2, pClippedWidth, NON_CODING_HEIGHT);

                } else {// Exon contains 5' UTR -- draw non-coding part
                    if (pCdStart > pStart) {
                        int pClippedStart = (int) Math.max(pStart, trackRectangle.getX());
                        int pClippedEnd = (int) Math.min(pCdStart, trackRectangle.getMaxX());
                        int pClippedWidth = pClippedEnd - pClippedStart;
                        drawExonRect(blockGraphics, exon, pClippedStart, curYOffset - NON_CODING_HEIGHT / 2, pClippedWidth, NON_CODING_HEIGHT);
                        pStart = pCdStart;
                    }
                    //  Exon contains 3' UTR  -- draw non-coding part
                    if (pCdEnd < pEnd) {
                        int pClippedStart = (int) Math.max(pCdEnd, trackRectangle.getX());
                        int pClippedEnd = (int) Math.min(pEnd, trackRectangle.getMaxX());
                        int pClippedWidth = pClippedEnd - pClippedStart;
                        drawExonRect(blockGraphics, exon, pClippedStart, curYOffset - NON_CODING_HEIGHT / 2, pClippedWidth, NON_CODING_HEIGHT);
                        pEnd = pCdEnd;
                    }

                    // At least part of this exon is coding.  Draw the coding part.
                    if ((exon.getCdStart() < exon.getEnd()) && (exon.getCdEnd() > exon.getStart())) {
                        int pClippedStart = (int) Math.max(pStart, trackRectangle.getX());
                        int pClippedEnd = (int) Math.min(pEnd, trackRectangle.getMaxX());
                        int pClippedWidth = Math.max(2, pClippedEnd - pClippedStart);
                        drawExonRect(blockGraphics, exon, pClippedStart, curYOffset - blockHeight / 2, pClippedWidth, blockHeight);
                    }
                }

                Graphics2D whiteArrowGraphics = context.getGraphic2DForColor(Color.white);
                drawStrandArrows(gene.getStrand(), pStart + ARROW_SPACING / 2, pEnd, curYOffset, 0, mode,
                        trackRectangle, whiteArrowGraphics);

                if (locationScale < 0.25) {
                    labelAminoAcids(pStart, fontGraphics, theOrigin, context, gene, locationScale,
                            curYOffset, trackRectangle, idx);
                }

            }
            exonNumberGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_DEFAULT);
            fontGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_DEFAULT);
        }
    }

    /**
     * Draw an individual rectangle representing an individual exon
     *
     * @param x
     * @param y
     * @param width
     * @param height
     */
    protected void drawExonRect(Graphics blockGraphics, Exon exon, int x, int y, int width, int height) {
        blockGraphics.fillRect(x, y, width, height);
    }

    /**
     * @param strand
     * @param startX
     * @param endX
     * @param startY
     * @param angle
     * @param mode
     * @param g2D
     */
    protected void drawStrandArrows(Strand strand,
                                    int startX,
                                    int endX,
                                    int startY,
                                    double angle,
                                    Track.DisplayMode mode,
                                    Rectangle trackRectangle,
                                    Graphics2D g2D) {


        //Don't draw arrows if we don't have a strand
        if (!strand.equals(Strand.POSITIVE) && !strand.equals(Strand.NEGATIVE)) {
            return;
        }

        // Don't draw strand arrows for very small regions
        // Limit drawing to visible region, we don't really know the viewport pEnd,
        if ((endX - startX)  < 6) {
            return;
        }

        Graphics2D g = (Graphics2D) g2D.create();


        // Draw the directional arrows on the feature
        int sz = mode == Track.DisplayMode.EXPANDED ? 3 : 2;
        sz = strand.equals(Strand.POSITIVE) ? -sz : sz;
        final int asz = Math.abs(sz);

        /*
         We draw arrows in a translated and rotated frame.
         This is to deal with alternative splice
         */

        int delta = startX < trackRectangle.x ? (int) ((trackRectangle.x - startX) % ARROW_SPACING) : 0;
        startX = Math.max(startX, trackRectangle.x);
        endX = Math.min(endX, trackRectangle.x + trackRectangle.width);

        int distance = endX - startX;

        g.translate(startX - delta, startY);
        g.rotate(-angle);
        double endXInFrame = distance / Math.cos(angle);

        for (int ii = ARROW_SPACING / 2; ii < endXInFrame; ii += ARROW_SPACING) {

            g.drawLine(ii, 0, ii + sz, +asz);
            g.drawLine(ii, 0, ii + sz, -asz);
        }

    }

    final private int drawFeatureName(IGVFeature feature, Track.DisplayMode mode, int pixelStart, int pixelEnd,
                                      int lastFeatureEndedAtPixelX, Graphics2D g2D,
                                      int textBaselineY) {

        String name = feature.getName();
        if (name == null) {
            return lastFeatureEndedAtPixelX;
        }

        FontMetrics fm = g2D.getFontMetrics();
        int fontSize = fm.getFont().getSize();
        int nameWidth = fm.stringWidth(name);
        int nameStart = (pixelStart + pixelEnd - nameWidth) / 2;

        Rectangle2D sb = fm.getStringBounds(name, g2D);


        if (nameStart > (lastFeatureEndedAtPixelX + fontSize) && sb.getWidth() < g2D.getClipBounds().getWidth()) {

            // g2D.clearRect(xString2, textBaselineY, (int) stringBounds.getWidth(), (int) stringBounds.getHeight());
            g2D.drawString(name, nameStart, textBaselineY);
            lastFeatureEndedAtPixelX = nameStart + nameWidth;
            drawnNames.add(name);

        }

        return lastFeatureEndedAtPixelX;
    }

    /**
     * Method description
     *
     * @param pStart
     * @param fontGraphics
     * @param theOrigin
     * @param context
     * @param gene
     * @param locationScale
     * @param yOffset
     * @param trackRectangle
     * @param idx  exon index
     */
    public void labelAminoAcids(int pStart, Graphics2D fontGraphics, double theOrigin,
                                RenderContext context, IGVFeature gene, double locationScale,
                                int yOffset, Rectangle trackRectangle, int idx) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        Exon exon = gene.getExons().get(idx);

        Exon prevExon = idx == 0 ? null : gene.getExons().get(idx-1);
        Exon nextExon = (idx+1) < gene.getExons().size() ?gene.getExons().get(idx+1) : null;

        AminoAcidSequence aaSequence = exon.getAminoAcidSequence(genome, prevExon, nextExon);

        if ((aaSequence != null) && aaSequence.hasNonNullSequence()) {
            Rectangle aaRect = new Rectangle(pStart, yOffset - blockHeight / 2, 1, blockHeight);

            int aaSeqStartPosition = aaSequence.getStartPosition();
            boolean odd =  exon.getAminoAcidNumber(exon.getCdStart()) % 2 == 1;

            for (AminoAcid acid : aaSequence.getSequence()) {
                if (acid != null) {

                    int start = Math.max(exon.getStart(), aaSeqStartPosition);
                    int end = Math.min(exon.getEnd(), aaSeqStartPosition + 3);
                    int px = getPixelFromChromosomeLocation(exon.getChr(), start, theOrigin, locationScale);
                    int px2 = getPixelFromChromosomeLocation(exon.getChr(), end, theOrigin, locationScale);

                    if ((px <= trackRectangle.getMaxX()) && (px2 >= trackRectangle.getX())) {

                        aaRect.x = px;
                        aaRect.width = px2 - px;


                        Graphics2D bgGraphics = context.getGraphic2DForColor(odd ? AA_COLOR_1 : AA_COLOR_2);
                        if (((acid.getSymbol() == 'M') && (((gene.getStrand() == Strand.POSITIVE) &&
                                (aaSeqStartPosition == exon.getCdStart())) || ((gene.getStrand() == Strand.NEGATIVE) &&
                                (aaSeqStartPosition == exon.getCdEnd() - 3))))) {
                            bgGraphics = context.getGraphic2DForColor(Color.green);
                        } else if (acid.getSymbol() == '*') {
                            bgGraphics = context.getGraphic2DForColor(Color.RED);
                        }

                        bgGraphics.fill(aaRect);

                        String tmp = new String(new char[]{acid.getSymbol()});
                        GraphicUtils.drawCenteredText(tmp, aaRect, fontGraphics);
                    }
                    odd = !odd;
                    aaSeqStartPosition += 3;

                }
            }

            if (aaSeqStartPosition < exon.getEnd()) {

                // The last codon is not complete (continues on next exon).
                // TODO -- get base from previous exon and compute aaSequence
                int cdEnd = Math.min(exon.getCdEnd(), exon.getEnd());
                aaRect.x = getPixelFromChromosomeLocation(exon.getChr(), aaSeqStartPosition, theOrigin, locationScale);
                aaRect.width = getPixelFromChromosomeLocation(exon.getChr(), cdEnd, theOrigin, locationScale) - aaRect.x;
                Graphics2D bgGraphics = context.getGraphic2DForColor(odd ? AA_COLOR_1 : AA_COLOR_2);
                odd = !odd;
                bgGraphics.fill(aaRect);
            }
        }
    }

    /**
     * Method description
     *
     * @return
     */
    public String getDisplayName() {
        return "Basic Feature";
    }

    protected Color getFeatureColor(IGVFeature feature, Track track) {
        // Set color used to draw the feature

        Color color = null;
        if (track.isItemRGB()) {
            color = feature.getColor();
        }

        if (color == null) {
            // TODO -- hack, generalize this
            if (track.getTrackType() == TrackType.CNV) {
                color = feature.getName().equals("gain") ? DULL_RED : DULL_BLUE;
            } else {
                // Only used if feature color is not set
                color = track.getColor();
            }

        }

        if (track.isUseScore()) {
            float min = getViewLimitMin(track);
            float max = getViewLimitMax(track);

            float score = feature.getScore();
            float alpha = Float.isNaN(score) ? 1 : getAlpha(min, max, feature.getScore());
            color = ColorUtilities.getCompositeColor(color, alpha);
        }

        return color;
    }

    float getViewLimitMin(Track track) {
        if (Float.isNaN(viewLimitMin)) {
            viewLimitMin = Float.isNaN(track.getViewLimitMin()) ? 0 : track.getViewLimitMin();
        }
        return viewLimitMin;
    }

    float getViewLimitMax(Track track) {
        if (Float.isNaN(viewLimitMax)) {
            viewLimitMax = Float.isNaN(track.getViewLimitMax()) ? 1000 : track.getViewLimitMax();
        }
        return viewLimitMax;
    }

    private float getAlpha(float minRange, float maxRange, float value) {
        float binWidth = (maxRange - minRange) / 9;
        int binNumber = (int) ((value - minRange) / binWidth);
        return Math.min(1.0f, 0.2f + (binNumber * 0.8f) / 9);
    }


    protected int getPixelFromChromosomeLocation(String chr, int chromosomeLocation, double origin,
                                                 double locationScale) {
        return (int) Math.round((chromosomeLocation - origin) / locationScale);
    }

    @Override
    public void reset() {
        exonMap.clear();
        drawnNames.clear();
    }
}
