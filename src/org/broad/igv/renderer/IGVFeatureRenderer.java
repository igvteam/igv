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
    static final Color AA_COLOR_1 = new Color(92, 92, 164);
    static final Color AA_COLOR_2 = new Color(12, 12, 120);
    static final Color DULL_BLUE = new Color(0, 0, 200);
    static final Color DULL_RED = new Color(200, 0, 0);
    static final int NON_CODING_HEIGHT = 8;

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
    private AlternativeSpliceGraph<Integer> exonGraph = new AlternativeSpliceGraph<Integer>();
    private Set<String> drawnNames = new HashSet<String>(100);

    protected boolean isGenotypeRenderer = false;


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
            fontGraphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

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

                int pixelEnd = (int) Math.round(Math.min(trackRectangleMaxX, virtualPixelEnd));
                int pixelStart = (int) Math.round(Math.max(trackRectangleX, virtualPixelStart));

                if (isGenotypeRenderer) {
                    if ((pixelEnd - pixelStart) < 3) {
                        double dx = 3.0 - (pixelEnd - pixelStart);
                        pixelStart -= dx / 2.0;
                        pixelEnd += dx / 2.0;
                    }
                }

                // Draw a maximum of "maxOcclusions" features on top of each other.
                if (pixelEnd <= lastPixelEnd) {
                    if (occludedCount >= maxOcclusions) {
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
                    pixelThickStart = (int) Math.max(trackRectangleX, (bf.getThickStart() - origin) / locScale);
                    pixelThickEnd = (int) Math.min(trackRectangleMaxX, (bf.getThickEnd() - origin) / locScale);
                    hasExons = bf.hasExons();
                }

                // Add directional arrows and exons, if there is room.
                int pixelYCenter = trackRectangle.y + NORMAL_STRAND_Y_OFFSET / 2;
                if (!hasExons) {
                    drawFeatureBlock(pixelStart, pixelEnd, pixelThickStart, pixelThickEnd, pixelYCenter, g2D);
                }


                if ((pixelEnd - pixelStart < 3) && hasExons) {
                    drawFeatureBounds(pixelStart, pixelEnd, pixelYCenter, g2D);
                } else {
                    if (hasExons) {
                        drawExons(feature, pixelYCenter, context, g2D, trackRectangle, displayMode,
                                alternateExonColor, track.getColor(), track.getAltColor());
                    } else {
                        Graphics2D arrowGraphics = hasExons
                                ? g2D
                                : context.getGraphic2DForColor(Color.WHITE);
                        drawStrandArrows(feature.getStrand(), pixelStart, pixelEnd, pixelYCenter, 0,
                                displayMode, arrowGraphics);
                    }

                }

                // Draw name , if there is room
                if (displayMode != Track.DisplayMode.SQUISHED) {
                    String name = feature.getName();
                    if (name != null) {
                        LineMetrics lineMetrics = font.getLineMetrics(name, g2D.getFontRenderContext());
                        int fontHeight = (int) Math.ceil(lineMetrics.getHeight());


                        // Draw feature name.  Center it over the feature extent,
                        // or if the feature extends beyond the track bounds over
                        // the track rectangle.
                        int nameStart = Math.max(0, pixelStart);
                        int nameEnd = Math.min(pixelEnd, (int) trackRectangle.getWidth());
                        int textBaselineY = pixelYCenter + blockHeight;

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
                    highlightGraphics.drawRect(pixelStart - 1, yStart, (pixelEnd - pixelStart + 2), blockHeight + 2);


                }


            }

            if (drawBoundary) {
                Graphics2D g2D = context.getGraphic2DForColor(Color.LIGHT_GRAY);
                g2D.drawLine((int) trackRectangleX, (int) trackRectangleMaxY - 1,
                        (int) trackRectangleMaxX, (int) trackRectangleMaxY - 1);
            }
        }
    }

    protected void renderGenotypeFeature(RenderContext context, IGVFeature feature, Rectangle trackRectangle, int pixelStart, int pixelEnd) {
        Color color = feature.getName().equals("HET") ? VariantRenderer.colorHet : VariantRenderer.colorHomVar;
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
    final private void drawFeatureBlock(int pixelStart, int pixelEnd, int pixelThickStart, int pixelThickEnd,
                                        int yOffset, Graphics2D g) {

        Graphics2D g2D = (Graphics2D) g.create();
        g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

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
        g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        float lineThickness = ((BasicStroke) g.getStroke()).getLineWidth();
        if (strand == null) {

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

        Graphics exonNumberGraphics = g2D.create();
        ((Graphics2D) exonNumberGraphics).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        exonNumberGraphics.setColor(Color.BLACK);
        exonNumberGraphics.setFont(FontManager.getFont(Font.BOLD, 8));

        // Now get the individual regions of the
        // feature are drawn here

        double theOrigin = context.getOrigin();
        double locationScale = context.getScale();

        Graphics2D fontGraphics = context.getGraphic2DForColor(Color.WHITE);

        boolean colorToggle = true;

        int lastX = Integer.MIN_VALUE;
        int lastY = Integer.MIN_VALUE;
        IExon lastExon = null;

        exonGraph.startFeature();
        for (Exon exon : gene.getExons()) {

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

            if (mode == Track.DisplayMode.ALTERNATIVE_SPLICE) {

                IExon eProx = Exon.getExonProxy(exon);
                drawConnectingLine = !exonGraph.containsEdge(lastExon, eProx);
                if (exonGraph.hasParameter(eProx)) {
                    curYOffset = exonGraph.getParameter(eProx);
                } else {
                    exonGraph.put(eProx, curYOffset);
                }
                lastExon = eProx;
            }

            int pStart = getPixelFromChromosomeLocation(exon.getChr(), exon.getStart(), theOrigin, locationScale);
            int pEnd = getPixelFromChromosomeLocation(exon.getChr(), exon.getEnd(), theOrigin, locationScale);


            Graphics2D arrowGraphics = context.getGraphic2DForColor(Color.blue);
            if (drawConnectingLine && lastX > Integer.MIN_VALUE && lastY > Integer.MIN_VALUE) {
                drawConnectingLine(lastX, lastY, pStart, curYOffset, exon.getStrand(), blockGraphics);
                double angle = Math.atan(-(curYOffset - lastY) / ((pStart - lastX) + 1e-12));
                drawStrandArrows(gene.getStrand(), lastX, pStart, lastY, angle, mode, arrowGraphics);
            }
            lastX = pEnd;
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
                if (exon.isUTR()) {
                    int pClippedStart = (int) Math.max(pStart, trackRectangle.getX());
                    int pClippedEnd = (int) Math.min(pEnd, trackRectangle.getMaxX());
                    int pClippedWidth = pClippedEnd - pClippedStart;
                    blockGraphics.fillRect(pClippedStart, curYOffset - NON_CODING_HEIGHT / 2, pClippedWidth, NON_CODING_HEIGHT);

                } else {// Exon contains 5' UTR -- draw non-coding part
                    if (pCdStart > pStart) {
                        int pClippedStart = (int) Math.max(pStart, trackRectangle.getX());
                        int pClippedEnd = (int) Math.min(pCdStart, trackRectangle.getMaxX());
                        int pClippedWidth = pClippedEnd - pClippedStart;
                        blockGraphics.fillRect(pClippedStart, curYOffset - NON_CODING_HEIGHT / 2, pClippedWidth, NON_CODING_HEIGHT);
                        pStart = pCdStart;
                    }
                    //  Exon contains 3' UTR  -- draw non-coding part
                    if (pCdEnd < pEnd) {
                        int pClippedStart = (int) Math.max(pCdEnd, trackRectangle.getX());
                        int pClippedEnd = (int) Math.min(pEnd, trackRectangle.getMaxX());
                        int pClippedWidth = pClippedEnd - pClippedStart;
                        blockGraphics.fillRect(pClippedStart, curYOffset - NON_CODING_HEIGHT / 2, pClippedWidth, NON_CODING_HEIGHT);
                        pEnd = pCdEnd;
                    }

                    // At least part of this exon is coding.  Draw the coding part.
                    if ((exon.getCdStart() < exon.getEnd()) && (exon.getCdEnd() > exon.getStart())) {
                        int pClippedStart = (int) Math.max(pStart, trackRectangle.getX());
                        int pClippedEnd = (int) Math.min(pEnd, trackRectangle.getMaxX());
                        int pClippedWidth = Math.max(2, pClippedEnd - pClippedStart);
                        blockGraphics.fillRect(pClippedStart, curYOffset - blockHeight / 2, pClippedWidth, blockHeight);

                    }
                }

                Graphics2D whiteArrowGraphics = context.getGraphic2DForColor(Color.white);
                drawStrandArrows(gene.getStrand(), pStart + ARROW_SPACING / 2, pEnd, curYOffset, 0, mode,
                        whiteArrowGraphics);

                if (locationScale < 0.25) {
                    labelAminoAcids(pStart, fontGraphics, theOrigin, context, gene, locationScale,
                            curYOffset, exon, trackRectangle);
                }

                if (FrameManager.isExomeMode()) {
                    //Draw exon edge bound
                    int halfHeight = blockHeight / 2;
                    edgeGraphics.drawLine(pStart, curYOffset - halfHeight, pStart, curYOffset + halfHeight - 1);
                    edgeGraphics.drawLine(pEnd, curYOffset - halfHeight, pEnd, curYOffset + halfHeight - 1);
                }

            }

        }

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
    protected void drawStrandArrows(Strand strand, int startX, int endX, int startY, double angle, Track.DisplayMode mode,
                                    Graphics2D g2D) {

        // Don't draw strand arrows for very small regions

        int distance = endX - startX;
        if ((distance < 6)) {
            return;
        }

        Graphics2D g = (Graphics2D) g2D.create();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        // Limit drawing to visible region, we don't really know the viewport pEnd,
        int vStart = 0;
        int vEnd = 10000;

        // Draw the directional arrows on the feature
        int sz = mode == Track.DisplayMode.EXPANDED ? 3 : 2;
        sz = strand.equals(Strand.POSITIVE) ? -sz : sz;

        final int asz = Math.abs(sz);

        /*
         We draw arrows in a translated and rotated frame.
         This is to deal with alternative splice
         */
        g.translate(startX, startY);
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
        if (name == null || (drawnNames.contains(name) && mode == Track.DisplayMode.ALTERNATIVE_SPLICE)) {
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
     * @param exon
     * @param trackRectangle
     */
    public void labelAminoAcids(int pStart, Graphics2D fontGraphics, double theOrigin,
                                RenderContext context, IGVFeature gene, double locationScale,
                                int yOffset, Exon exon, Rectangle trackRectangle) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        AminoAcidSequence aaSequence = exon.getAminoAcidSequence(genome);
        if ((aaSequence != null) && aaSequence.hasNonNullSequence()) {
            Rectangle aaRect = new Rectangle(pStart, yOffset - blockHeight / 2, 1, blockHeight);

            int aaSeqStartPosition = aaSequence.getStartPosition();
            int firstFullAcidIndex = (int) Math.floor((aaSeqStartPosition - exon.getReadingShift()) / 3);
            //calculated oddness or evenness of first amino acid. This is also done independently in SequenceRenderer
            boolean odd = (firstFullAcidIndex % 2) == 1;

            if (aaSeqStartPosition > exon.getStart()) {

                // The codon for the first amino acid is split between this and the previous exon
                // TODO -- get base from previous exon and compute aaSequence.  For now skipping drawing
                // the AA label
                int cdStart = Math.max(exon.getCdStart(), exon.getStart());
                aaRect.x = getPixelFromChromosomeLocation(exon.getChr(), cdStart, theOrigin,
                        locationScale);
                aaRect.width = getPixelFromChromosomeLocation(exon.getChr(),
                        aaSeqStartPosition, theOrigin, locationScale) - aaRect.x;

                if (trackRectangle.intersects(aaRect)) {
                    //use opposite color from first AA color
                    Graphics2D bgGraphics = context.getGraphic2DForColor(!odd ? AA_COLOR_1 : AA_COLOR_2);
                    bgGraphics.fill(aaRect);
                }
            }


            for (AminoAcid acid : aaSequence.getSequence()) {
                if (acid != null) {

                    int px = getPixelFromChromosomeLocation(exon.getChr(), aaSeqStartPosition, theOrigin, locationScale);
                    int px2 = getPixelFromChromosomeLocation(exon.getChr(), aaSeqStartPosition + 3, theOrigin, locationScale);

                    if ((px <= trackRectangle.getMaxX()) && (px2 >= trackRectangle.getX())) {

                        // {

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
        exonGraph = new AlternativeSpliceGraph<Integer>();
        drawnNames.clear();
    }
}
