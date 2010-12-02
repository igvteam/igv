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
package org.broad.igv.renderer;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.FontManager;
import org.broad.igv.util.ColorUtilities;
import org.broad.tribble.Feature;

import java.awt.*;
import java.awt.font.LineMetrics;
import java.util.List;

/**
 * Renderer for splice junctions
 *
 * @author Enter your name here...
 * @version Enter version here..., 08/10/24
 */
public class SpliceJunctionRenderer extends IGVFeatureRenderer {

    private static Logger log = Logger.getLogger(SpliceJunctionRenderer.class);


    //color for drawing all arcs
    Color ARC_COLOR_POS = new Color(50, 50, 150, 140);
    Color ARC_COLOR_NEG = new Color(150, 50, 50, 140);


    //maximum depth that can be displayed
    protected int MAX_DEPTH = 50;


    /**
     * Note:  assumption is that featureList is sorted by pStart position.
     *
     * @param featureList
     * @param context
     * @param trackRectangle
     * @param track
     */
    public void renderFeatures(List<IGVFeature> featureList,
                               RenderContext context,
                               Rectangle trackRectangle,
                               FeatureTrack track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        // Clear values
        lastFeatureLineMaxY = 0;
        lastFeatureBoundsMaxY = 0;
        lastRegionMaxY = 0;

        // TODO -- use enum instead of string "Color"
        if ((featureList != null) && (featureList.size() > 0)) {

            // Create a graphics object to draw font names.  Graphics are not cached
            // by font, only by color, so its neccessary to create a new one to prevent
            // affecting other tracks.
            font = FontManager.getScalableFont(track.getFontSize());
            Graphics2D fontGraphics = (Graphics2D) context.getGraphic2DForColor(Color.BLACK).create();
            fontGraphics.setFont(font);

            // Track coordinates
            double trackRectangleX = trackRectangle.getX();
            double trackRectangleMaxX = trackRectangle.getMaxX();
            double trackRectangleY = trackRectangle.getY();
            double trackRectangleMaxY = trackRectangle.getMaxY();

            // Draw the lines that represent the bounds of
            // a feature's region
            // TODO -- bugs in "Line Placement" style -- hardocde to fishbone


            int lastFeatureEndedAtPixelX = -9999;
            int lastPixelEnd = -1;
            int occludedCount = 0;
            int maxOcclusions = 2;

            IGVFeature featureArray[] = featureList.toArray(new IGVFeature[featureList.size()]);
            for (int i = 0; i < featureArray.length; i++) {

                IGVFeature feature = featureArray[i];


                // Get the pStart and pEnd of the entire feature  at extreme zoom levels the
                // virtual pixel value can be too large for an int, so the computation is
                // done in double precision and cast to an int only when its confirmed its
                // within the field of view.


                double virtualPixelStart = Math.round((feature.getStart() - origin) / locScale);
                double virtualPixelEnd = Math.round((feature.getEnd() - origin) / locScale);

                boolean hasRegions = ((feature.getExons() != null) && (feature.getExons().size() > 0));

                // If the any part of the feature fits in the
                // Track rectangle draw it
                if ((virtualPixelEnd >= trackRectangleX) && (virtualPixelStart <= trackRectangleMaxX)) {

                    //
                    int pixelEnd = (int) virtualPixelEnd;//(int) Math.min(trackRectangleMaxX, virtualPixelEnd);
                    int pixelStart = (int) virtualPixelStart;//(int) Math.max(trackRectangleX, virtualPixelStart);


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

                    Color color = ARC_COLOR_POS;
                    if (feature.getStrand() != null && feature.getStrand().equals(Strand.NEGATIVE))
                        color = ARC_COLOR_NEG;


                    Graphics2D g2D = context.getGraphic2DForColor(color);

                    // Get the feature's direction
                    Strand strand = feature.getStrand();

                    // Draw block
                    float score = 1;
                    if (feature instanceof BasicFeature) {
                        BasicFeature bf = (BasicFeature) feature;
                        score = bf.getScore();
                    }

                    int pixelScore = Math.max(1,(int) ((Math.min(score, MAX_DEPTH) / MAX_DEPTH) * trackRectangle.getHeight()));

                    drawFeatureBlock(pixelStart, pixelEnd, pixelScore, trackRectangle, strand, hasRegions, g2D);

                    // Determine the y offset of features based on strand type

                    // If the width is < 3 there isn't room to draw the
                    // feature, or orientation.  If the feature has any exons
                    // at all indicate by filling a small rect

                    int pixelYCenter = trackRectangle.y + NORMAL_STRAND_Y_OFFSET / 2;



                    String name = feature.getName();
                    if (name != null) {
                        LineMetrics lineMetrics = font.getLineMetrics(name, g2D.getFontRenderContext());

                    }

                    // If this is the highlight feature highlight it
                    if (getHighlightFeature() == feature) {
                        int yStart = pixelYCenter - BLOCK_HEIGHT / 2 - 1;
                        Graphics2D highlightGraphics = context.getGraphic2DForColor(Color.cyan);
                        highlightGraphics.drawRect(pixelStart - 1, yStart, (pixelEnd - pixelStart + 2), BLOCK_HEIGHT + 2);


                    }

                }

            }

            if (drawBoundary) {
                Graphics2D g2D = context.getGraphic2DForColor(Color.LIGHT_GRAY);
                g2D.drawLine((int) trackRectangleX, (int) trackRectangleMaxY - 1,
                        (int) trackRectangleMaxX, (int) trackRectangleMaxY - 1);
            }
        }
    }

    /**
     *
     * @param pixelStart
     * @param pixelEnd
     * @param pixelScore
     * @param trackRectangle
     * @param strand
     * @param hasRegions
     * @param g
     */
    final private void drawFeatureBlock(int pixelStart, int pixelEnd, int pixelScore,
                                        Rectangle trackRectangle,
                                        Strand strand, boolean hasRegions, Graphics2D g) {

        Graphics2D g2D = (Graphics2D) g.create();

        int yOffset = trackRectangle.y + NORMAL_STRAND_Y_OFFSET / 2;

        int pixelWidth = Math.max(1,pixelEnd - pixelStart);
        for (int i=0; i<pixelScore; i++)
        {
            g2D.drawArc(pixelStart, yOffset + i, pixelWidth, (int) trackRectangle.getHeight() - 2*i, 0, 180);
        }

        lastFeatureLineMaxY = yOffset + pixelScore;

        g2D.dispose();
    }


}
