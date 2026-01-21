 package org.igv.renderer;


 import htsjdk.tribble.Feature;
 import org.igv.Globals;
 import org.igv.logging.*;
 import org.igv.feature.IGVFeature;
 import org.igv.feature.SpliceJunctionFeature;
 import org.igv.feature.Strand;
 import org.igv.prefs.Constants;
 import org.igv.prefs.PreferencesManager;
 import org.igv.track.FeatureTrack;
 import org.igv.track.RenderContext;
 import org.igv.track.Track;
 import org.igv.ui.color.ColorUtilities;

 import java.awt.*;
 import java.awt.geom.GeneralPath;
 import java.util.List;

 /**
  * Renderer for splice junctions. Draws a filled-in arc for each junction, with the width of the
  * arc representing the depth of coverage. If coverage information is present for the flanking
  * regions, draws that, too; otherwise indicates flanking regions with rectangles
  *
  * @author dhmay
  */
 public class SpliceJunctionRenderer extends IGVFeatureRenderer {

     private static Logger log = LogManager.getLogger(SpliceJunctionRenderer.class);

     //color for drawing all arcs
     private static Color ARC_COLOR_NEG = new Color(50, 50, 150, 140); //transparent dull blue
     private static Color ARC_COLOR_POS = new Color(150, 50, 50, 140); //transparent dull red
     private static Color ARC_COLOR_HIGHLIGHT_NEG = new Color(90, 90, 255, 255); //opaque, brighter blue
     private Color ARC_COLOR_HIGHLIGHT_POS = new Color(255, 90, 90, 255); //opaque, brighter red
     private static Color COLOR_CENTERLINE = new Color(0, 0, 0, 100);

     private int maxDepth = 50;

     public SpliceJunctionRenderer() {
         if (Globals.isDarkMode()) {
             // Use brighter colors in dark mode
             ARC_COLOR_NEG = ColorUtilities.modifyAlpha(Color.CYAN, 140);
             ARC_COLOR_POS = ColorUtilities.modifyAlpha(Color.RED, 140);
             ARC_COLOR_HIGHLIGHT_NEG = Color.BLUE;
             ARC_COLOR_HIGHLIGHT_POS = Color.RED;
             COLOR_CENTERLINE = new Color(255, 255, 255, 100);
         }
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


         // TODO -- use enum instead of string "Color"
         if ((featureList != null) && !featureList.isEmpty()) {

             Graphics2D g2D = null;

             try {
                 g2D = (Graphics2D) context.getGraphics().create();

                 double origin = context.getOrigin();
                 double locScale = context.getScale();
                 double end = origin + locScale * trackRectangle.width;

                 Feature selectedFeature = ((FeatureTrack) track).getSelectedFeature();
                 boolean shouldShowFlankingRegions = PreferencesManager.getPreferences().getAsBoolean(Constants.SAM_SHOW_JUNCTION_FLANKINGREGIONS);

                 for (IGVFeature feature : featureList) {
                     SpliceJunctionFeature junctionFeature = (SpliceJunctionFeature) feature;

                     // If any part of the feature fits in the track rectangle draw it
                     if (junctionFeature.getEnd() > origin && junctionFeature.getStart() < end) {
                         boolean shouldHighlight = junctionFeature.isSameJunction(selectedFeature);
                         drawFeature(junctionFeature, shouldShowFlankingRegions, shouldHighlight, origin, locScale, trackRectangle, track, g2D);
                     }
                 }

                 //draw a central horizontal line
                 g2D.setColor(COLOR_CENTERLINE);
                 g2D.drawLine(trackRectangle.x, (int) trackRectangle.getCenterY(),
                         trackRectangle.x + trackRectangle.width, (int) trackRectangle.getCenterY());
             } finally {
                 g2D.dispose();
             }

         }
     }

     /**
      * Draw a filled arc representing a single feature. The thickness and height of the arc are proportional to the
      * depth of coverage.  Some of this gets a bit arcane -- the result of lots of visual tweaking.
      *
      * @param junctionFeature
      * @param highlight
      * @param trackRectangle
      * @param g2D
      */
     protected void drawFeature(SpliceJunctionFeature junctionFeature, boolean shouldShowFlankingRegions,
                                boolean highlight, double origin, double locScale, Rectangle trackRectangle, Track track, Graphics2D g2D) {

         int flankingStart = junctionFeature.getStart();
         int flankingEnd = junctionFeature.getEnd();

         int junctionStart = junctionFeature.getJunctionStart();
         int junctionEnd = junctionFeature.getJunctionEnd();

         int pixelFeatureStart = (int) Math.round((flankingStart - origin) / locScale);
         int pixelFeatureEnd = (int) Math.round((flankingEnd - origin) / locScale);

         int pixelJunctionStart = (int) Math.round((junctionStart - origin) / locScale);
         int pixelJunctionEnd = (int) Math.round((junctionEnd - origin) / locScale);

         Strand strand = junctionFeature.getStrand();
         float depth = junctionFeature.getJunctionDepth();
         Color featureColor = junctionFeature.getColor();


         boolean isPositiveStrand = true;
         // Get the feature's direction, color appropriately
         if (strand != null && strand.equals(Strand.NEGATIVE))
             isPositiveStrand = false;

         /*
            Choose the feature color:
            1. Check if the user specified a positive / negative track color
            2. Check if the feature has its own color
            3. Use the default
            We modify the alpha value of the chosen color to indicate if it is highlighted
          */
         final Color color;
         Color trackColor = isPositiveStrand ? track.getColor() : track.getAltColor();
         if(trackColor == null) {
             trackColor = track.getDefaultColor();
         }
         if (trackColor != null) {
             color = adjustAlpha(trackColor, highlight);
         } else if (featureColor != null) {
             color = adjustAlpha(featureColor, highlight);
         } else {
             if (isPositiveStrand) {
                 color = highlight ? ARC_COLOR_HIGHLIGHT_POS : ARC_COLOR_POS;
             } else {
                 color = highlight ? ARC_COLOR_HIGHLIGHT_NEG : ARC_COLOR_NEG;
             }
         }

         g2D.setColor(color);

         //Height of top of an arc of maximum depth
         int maxPossibleArcHeight = (trackRectangle.height - 1) / 2;

         if (shouldShowFlankingRegions) {
             if (junctionFeature.hasFlankingRegionDepthArrays()) {
                 //draw a wigglegram of the splice junction flanking region depth of coverage

                 int startFlankingRegionPixelLength = pixelJunctionStart - pixelFeatureStart;
                 int endFlankingRegionPixelLength = pixelFeatureEnd - pixelJunctionEnd;

                 drawFlankingRegion(g2D, pixelFeatureStart, startFlankingRegionPixelLength,
                         junctionFeature.getStartFlankingRegionDepthArray(), maxPossibleArcHeight,
                         trackRectangle, isPositiveStrand);
                 drawFlankingRegion(g2D, pixelJunctionEnd + 1, endFlankingRegionPixelLength,
                         junctionFeature.getEndFlankingRegionDepthArray(), maxPossibleArcHeight,
                         trackRectangle, isPositiveStrand);
             } else {
                 //Draw rectangles indicating the overlap on each side of the junction
                 int overlapRectHeight = 3;
                 int overlapRectTopX = (int) trackRectangle.getCenterY() + (isPositiveStrand ? -2 : 0);
                 if (pixelFeatureStart < pixelJunctionStart) {
                     g2D.fillRect(pixelFeatureStart, overlapRectTopX,
                             pixelJunctionStart - pixelFeatureStart, overlapRectHeight);
                 }
                 if (pixelJunctionEnd < pixelFeatureEnd) {
                     g2D.fillRect(pixelJunctionEnd, overlapRectTopX,
                             pixelFeatureEnd - pixelJunctionEnd, overlapRectHeight);
                 }
             }
         }

         //Create a path describing the arc, using Bezier curves. The Bezier control points for the top and
         //bottom arcs are based on the boundary points of the rectangles containing the arcs

         //proportion of the maximum arc height used by a minimum-height arc
         double minArcHeightProportion = 0.33;

         int innerArcHeight = (int) (maxPossibleArcHeight * minArcHeightProportion);
         float depthProportionOfMax = Math.min(1, depth / maxDepth);
         int arcWidth = Math.max(1, (int) ((1 - minArcHeightProportion) * maxPossibleArcHeight * depthProportionOfMax));
         int outerArcHeight = innerArcHeight + arcWidth;


         //Height of bottom of the arc
         int arcBeginY = (int) trackRectangle.getCenterY() + (isPositiveStrand ? -1 : 1);
         int outerArcPeakY = isPositiveStrand ? arcBeginY - outerArcHeight : arcBeginY + outerArcHeight;
         int innerArcPeakY = isPositiveStrand ? arcBeginY - innerArcHeight : arcBeginY + innerArcHeight;

         //dhmay: I don't really understand Bezier curves.  For some reason I have to put the Bezier control
         //points farther up or down than I want the arcs to extend.  This multiplier seems about right
         int outerBezierY = arcBeginY + (int) (1.3 * (outerArcPeakY - arcBeginY));
         int innerBezierY = arcBeginY + (int) (1.3 * (innerArcPeakY - arcBeginY));

         //Putting the Bezier control points slightly off to the sides of the arc
         int bezierXPad = Math.max(1, (pixelJunctionEnd - pixelJunctionStart) / 30);

         GeneralPath arcPath = new GeneralPath();
         arcPath.moveTo(pixelJunctionStart, arcBeginY);
         arcPath.curveTo(pixelJunctionStart - bezierXPad, outerBezierY, //Bezier 1
                 pixelJunctionEnd + bezierXPad, outerBezierY,         //Bezier 2
                 pixelJunctionEnd, arcBeginY);        //Arc end
         arcPath.curveTo(pixelJunctionEnd + bezierXPad, innerBezierY, //Bezier 1
                 pixelJunctionStart - bezierXPad, innerBezierY,         //Bezier 2
                 pixelJunctionStart, arcBeginY);        //Arc end

         g2D.setClip(new Rectangle(pixelJunctionStart, trackRectangle.y, pixelJunctionEnd - pixelJunctionStart, trackRectangle.height));

         //Draw the arc, to ensure outline is drawn completely (fill won't do it, necessarily). This will also
         //give the arc a darker outline
         g2D.draw(arcPath);
         //Fill the arc
         g2D.fill(arcPath);

         g2D.setClip(null);
         //g2D.drawLine(pixelJunctionStart, trackRectangle.y, pixelJunctionStart, trackRectangle.y + trackRectangle.height);
         //g2D.drawLine(pixelJunctionEnd, trackRectangle.y, pixelJunctionEnd, trackRectangle.y + trackRectangle.height);

     }

     /**
      * @return a variant of the input color with a different alpha depending on if it is a highlight color or not
      */
     private static Color adjustAlpha(final Color featureColor, final boolean highlight) {
         Color color;
         int r = featureColor.getRed();
         int g = featureColor.getGreen();
         int b = featureColor.getBlue();
         int alpha = highlight ? 255 : 140;
         color = new Color(r, g, b, alpha);
         return color;
     }


     /**
      * Draw depth of coverage for the starting or ending flanking region
      *
      * @param g2D
      * @param pixelStart
      * @param pixelLength
      * @param regionDepthArray
      * @param maxPossibleArcHeight
      * @param trackRectangle
      * @param isPositiveStrand
      */
     protected void drawFlankingRegion(Graphics g2D, int pixelStart, int pixelLength, int[] regionDepthArray,
                                       int maxPossibleArcHeight, Rectangle trackRectangle, boolean isPositiveStrand) {
         for (int i = 0; i < pixelLength; i++) {
             float arrayIndicesPerPixel = (float) regionDepthArray.length /
                     (float) pixelLength;
             int flankingRegionArrayPixelMinIndex = (int) (i * arrayIndicesPerPixel);
             int flankingRegionArrayPixelMaxIndex = (int) ((i + 1) * arrayIndicesPerPixel);
             flankingRegionArrayPixelMinIndex =
                     Math.max(0, Math.min(flankingRegionArrayPixelMinIndex, regionDepthArray.length - 1));
             flankingRegionArrayPixelMaxIndex =
                     Math.max(0, Math.min(flankingRegionArrayPixelMaxIndex, regionDepthArray.length - 1));

             int meanDepthThisPixel = 0;
             for (int j = flankingRegionArrayPixelMinIndex; j <= flankingRegionArrayPixelMaxIndex; j++)
                 meanDepthThisPixel += regionDepthArray[j];
             meanDepthThisPixel /= (flankingRegionArrayPixelMaxIndex - flankingRegionArrayPixelMinIndex + 1);
             meanDepthThisPixel = Math.min(maxDepth, meanDepthThisPixel);
             int pixelHeight = Math.max(maxPossibleArcHeight * meanDepthThisPixel / maxDepth, 2);
             g2D.fillRect(pixelStart + i,
                     (int) trackRectangle.getCenterY() + (isPositiveStrand ? -pixelHeight : 0),
                     1, pixelHeight);
         }
     }


     public int getMaxDepth() {
         return maxDepth;
     }

     public void setMaxDepth(int maxDepth) {
         this.maxDepth = maxDepth;
     }
 }
