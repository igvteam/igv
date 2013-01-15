/**
 * Copyright (c) 2010-2011 by Fred Hutchinson Cancer Research Center.  All Rights Reserved.

 * This software is licensed under the terms of the GNU Lesser General
 * Public License (LGPL), Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.

 * THE SOFTWARE IS PROVIDED "AS IS." FRED HUTCHINSON CANCER RESEARCH CENTER MAKES NO
 * REPRESENTATIONS OR WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED,
 * INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS,
 * WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL FRED HUTCHINSON CANCER RESEARCH
 * CENTER OR ITS TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR
 * ANY DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR
 * CONSEQUENTIAL DAMAGES, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS,
 * REGARDLESS OF  WHETHER FRED HUTCHINSON CANCER RESEARCH CENTER SHALL BE ADVISED,
 * SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */
package org.broad.igv.renderer;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.tribble.Feature;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.util.List;

/**
 * Renderer for splice junctions. Draws a filled-in arc for each junction, with the width of the
 * arc representing the depth of coverage. If coverage information is present for the flanking
 * regions, draws that, too; otherwise indicates flanking regions with rectangles
 *
 * @author dhmay
 */
public class SashimiJunctionRenderer extends IGVFeatureRenderer {

    private static Logger log = Logger.getLogger(SashimiJunctionRenderer.class);

    //color for drawing all arcs
    Color ARC_COLOR_NEG = new Color(50, 50, 150, 140); //transparent dull blue
    Color ARC_COLOR_POS = new Color(150, 50, 50, 140); //transparent dull red

    Color ARC_COLOR_HIGHLIGHT_NEG = new Color(90, 90, 255, 255); //opaque, brighter blue
    Color ARC_COLOR_HIGHLIGHT_POS = new Color(255, 90, 90, 255); //opaque, brighter red


    //central horizontal line color
    Color COLOR_CENTERLINE = new Color(0, 0, 0, 100);

    //maximum depth that can be displayed, due to track height limitations. Junctions with
    //this depth and deeper will all look the same
    protected int MAX_DEPTH = 50;

    private ShapeType shapeType = ShapeType.ELLIPSE;

    public enum ShapeType{
        CIRCLE,
        ELLIPSE,
        TEXT
    }

    public void setShapeType(ShapeType shapeType){
        this.shapeType = shapeType;
    }

    public ShapeType getShapeType() {
        return shapeType;
    }

    /**
     * Note:  assumption is that featureList is sorted by pStart position.
     *
     * @param featureList
     * @param context
     * @param trackRectangle
     * @param track
     */
    @Override
    public void render(List<IGVFeature> featureList,
                       RenderContext context,
                       Rectangle trackRectangle,
                       Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();

        // TODO -- use enum instead of string "Color"
        if ((featureList != null) && !featureList.isEmpty()) {

            // Create a graphics object to draw font names.  Graphics are not cached
            // by font, only by color, so its neccessary to create a new one to prevent
            // affecting other tracks.
            Font font = FontManager.getFont(track.getFontSize());
            Graphics2D fontGraphics = (Graphics2D) context.getGraphic2DForColor(Color.BLACK).create();
            fontGraphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            fontGraphics.setFont(font);

            //determine whether to show flanking regions
            PreferenceManager prefs = PreferenceManager.getInstance();
            boolean shouldShowFlankingRegions = prefs.getAsBoolean(
            PreferenceManager.SAM_SHOW_JUNCTION_FLANKINGREGIONS);

            // Track coordinates
            double trackRectangleX = trackRectangle.getX();
            double trackRectangleMaxX = trackRectangle.getMaxX();

            // Draw the lines that represent the bounds of
            // a feature's region
            // TODO -- bugs in "Line Placement" style -- hardocde to fishbone

            Feature selectedFeature = ((FeatureTrack) track).getSelectedExon();

            boolean drawAbove = true;
            for (IGVFeature feature : featureList) {
                SpliceJunctionFeature junctionFeature = (SpliceJunctionFeature) feature;

                // Get the pStart and pEnd of the entire feature.  at extreme zoom levels the
                // virtual pixel value can be too large for an int, so the computation is
                // done in double precision and cast to an int only when its confirmed its
                // within the field of view.
                int flankingStart = junctionFeature.getStart();
                int flankingEnd = junctionFeature.getEnd();

                int junctionStart = junctionFeature.getJunctionStart();
                int junctionEnd = junctionFeature.getJunctionEnd();

                //Only show arcs for the selected feature, if applicable
                if (selectedFeature != null) {
                    if((junctionStart >= selectedFeature.getStart() && junctionStart <= selectedFeature.getEnd())
                            || (junctionEnd >= selectedFeature.getStart() && junctionEnd <= selectedFeature.getEnd())){

                    }else{
                        continue;
                    }
                }

                double virtualPixelStart = Math.round((flankingStart - origin) / locScale);
                double virtualPixelEnd = Math.round((flankingEnd - origin) / locScale);

                double virtualPixelJunctionStart = Math.round((junctionStart - origin) / locScale);
                double virtualPixelJunctionEnd = Math.round((junctionEnd - origin) / locScale);

                // If the any part of the feature fits in the
                // Track rectangle draw it
                if ((virtualPixelEnd >= trackRectangleX) && (virtualPixelStart <= trackRectangleMaxX)) {

                    //
                    int displayPixelEnd = (int) Math.min(trackRectangleMaxX, virtualPixelEnd);
                    int displayPixelStart = (int) Math.max(trackRectangleX, virtualPixelStart);

                    int depth = junctionFeature.getJunctionDepth();
                    Color color = feature.getColor();

                    drawFeature((int) virtualPixelStart, (int) virtualPixelEnd,
                            (int) virtualPixelJunctionStart, (int) virtualPixelJunctionEnd, depth,
                            trackRectangle, context, drawAbove, junctionFeature, color,
                            shouldShowFlankingRegions);
                    drawAbove = !drawAbove;
                }
            }

            //draw a central horizontal line
            Graphics2D g2D = context.getGraphic2DForColor(COLOR_CENTERLINE);
            g2D.drawLine((int) trackRectangleX, (int) trackRectangle.getCenterY(),
                    (int) trackRectangleMaxX, (int) trackRectangle.getCenterY());

        }
    }

    /**
     * Draw a filled arc representing a single feature. The thickness and height of the arc are proportional to the
     * depth of coverage.  Some of this gets a bit arcane -- the result of lots of visual tweaking.
     *
     * @param pixelFeatureStart  the starting position of the feature, whether on-screen or not
     * @param pixelFeatureEnd    the ending position of the feature, whether on-screen or not
     * @param pixelJunctionStart the starting position of the junction, whether on-screen or not
     * @param pixelJunctionEnd   the ending position of the junction, whether on-screen or not
     * @param depth              coverage depth
     * @param trackRectangle
     * @param context
     * @param drawAbove Whether to draw arc above or below midline
     * @param junctionFeature
     * @param featureColor       the color specified for this feature.  May be null.
     */
    protected void drawFeature(int pixelFeatureStart, int pixelFeatureEnd,
                               int pixelJunctionStart, int pixelJunctionEnd, int depth,
                               Rectangle trackRectangle, RenderContext context, boolean drawAbove,
                               SpliceJunctionFeature junctionFeature, Color featureColor,
                               boolean shouldShowFlankingRegions) {

        //If the feature color is specified, use it, except that we set our own alpha depending on whether
        //the feature is highlighted.  Otherwise default based on strand and highlight.
        Color color;
        if (featureColor != null) {
            int r = featureColor.getRed();
            int g = featureColor.getGreen();
            int b = featureColor.getBlue();
            int alpha = 140;
            color = new Color(r, g, b, alpha);
        } else {
            color = drawAbove ? ARC_COLOR_POS : ARC_COLOR_NEG;
        }

        Graphics2D g2D = context.getGraphic2DForColor(color);

        //Height of top of an arc of maximum depth
        int maxPossibleArcHeight = (trackRectangle.height - 1) / 2;

        //proportion of the maximum arc height used by a minimum-height arc
        double minArcHeightProportion = 0.1;

        float depthProportionOfMax = Math.min(1, depth / MAX_DEPTH);
        int arcHeight = Math.max(5, (int) ((1 - minArcHeightProportion) * maxPossibleArcHeight * depthProportionOfMax));

        //We adjust up or down depending on the strand
        int yStrandModifier = drawAbove ? -1 : 1;

        int arcBeginY = (int) trackRectangle.getCenterY() + yStrandModifier;


        //We use corners of a square as control points because why not
        //The control point is never actually reached
        int arcControlPeakY = arcBeginY + yStrandModifier * arcHeight;

        GeneralPath arcPath = new GeneralPath();
        arcPath.moveTo(pixelJunctionStart, arcBeginY);
        arcPath.curveTo(pixelJunctionStart, arcControlPeakY,
                        pixelJunctionEnd, arcControlPeakY,
                pixelJunctionEnd, arcBeginY);

        g2D.draw(arcPath);

        float midX = ((float) pixelJunctionStart + (float) pixelJunctionEnd) / 2;

//        //TODO Format number
//        Graphics2D strGraphics = context.getGraphic2DForColor(Color.black);
//        strGraphics.drawString("" + depth, midX, arcPeakY);

        double actArcPeakY = arcBeginY + yStrandModifier * Math.pow(0.5, 3) * (6) * arcHeight;

        //Draw shape to indicate depth
        float maxPossibleShapeHeight = maxPossibleArcHeight / 2;

        Shape shape = null;
        switch(shapeType){
            case CIRCLE:
                shape = createDepthCircle(maxPossibleShapeHeight, depthProportionOfMax, midX, actArcPeakY);
                break;
            case ELLIPSE:
                shape = createDepthEllipse(maxPossibleShapeHeight, depthProportionOfMax, midX, actArcPeakY);
                break;
            case TEXT:
                g2D.drawString("" + depth, midX, arcControlPeakY);
        }

        if(shape != null){
            g2D.draw(shape);
            g2D.fill(shape);
        }
    }

    private Shape createDepthEllipse(double maxPossibleShapeHeight, double depthProportionOfMax, double arcMidX, double actArcPeakY){
        double w = 5f;
        double x = arcMidX - w/2;

        double h = maxPossibleShapeHeight * depthProportionOfMax;

        //The ellipse is always specified from the top left corner
        double y = actArcPeakY - h / 2;

        return new Ellipse2D.Double(x, y, w, h);
    }

    private Shape createDepthCircle(double maxPossibleShapeHeight, double depthProportionOfMax, double arcMidX, double actArcPeakY){

        double h = maxPossibleShapeHeight * Math.sqrt(depthProportionOfMax);
        double w = h;
        double x = arcMidX - w/2;

        //The ellipse is always specified from the top left corner
        double y = actArcPeakY - h / 2;

        return new Ellipse2D.Double(x, y, w, h);
    }


}
