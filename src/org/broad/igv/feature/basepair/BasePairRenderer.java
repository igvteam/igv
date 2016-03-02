// previously ArcRenderer

package org.broad.igv.feature.basepair;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.awt.geom.GeneralPath;

// TODO: add vertical scaling option so arcs fit on track, or clip arcs outside of track rect?

public class BasePairRenderer {

    private static Logger log = Logger.getLogger(BasePairRenderer.class);

    Color ARC_COLOR_A = new Color(50, 50, 150, 140); //transparent dull blue
    Color ARC_COLOR_B = new Color(150, 50, 50, 140); //transparent dull red
    Color ARC_COLOR_C = new Color(50, 0, 50, 250);

    int dir = -1; // 1 for up, -1 for down

    // central horizontal line color
    Color COLOR_CENTERLINE = new Color(0, 0, 0, 100);

    public int getDirection(){
        return dir;
    }

    public void setDirection(int d){
        dir = d;
    }

    public void draw(BasePairData data, RenderContext context, Rectangle trackRectangle){

        double nucsPerPixel = context.getScale();
        double origin = context.getOrigin();
        //The location of the first base that is loaded
        int start = Math.max(0, (int) origin - 1);
        //The location of the last base that is loaded
        int end = (int) (origin + trackRectangle.width * nucsPerPixel) + 1;
        if (end <= start) return;

        // TODO: should make this a function (ugh no multiple value return from java functions)

        java.util.List<BasePairFeature> featureList = data.getFeatures(context.getChr());

        if (featureList != null) {

            for (BasePairFeature feature : featureList) {

                if(feature.startLeft > context.getEndLocation()) break;
                else if(feature.endRight < context.getOrigin()) continue;

                //System.out.println("Color: "+data.colors[i]);
                int arcCount = 0;
                // TODO: only render arcs whose interior overlaps the viewing area - i.e. if an arc starts left of the window and ends right of the window, should still render


                //System.out.println("In arcRenderer.draw(): ");
                //System.out.println("    track.width = " + trackRectangle.width);
                //System.out.println("    nucsPerPixel = " + nucsPerPixel);
                //System.out.println("    origin = " + origin);
                //System.out.println("    start = " + start);
                //System.out.println("    end = " + end);

                double startLeftPix = (feature.startLeft - origin) / nucsPerPixel;
                double startRightPix = (feature.startRight + 1.0 - origin) / nucsPerPixel;
                double endLeftPix = (feature.endLeft - origin) / nucsPerPixel;
                double endRightPix = (feature.endRight + 1.0 - origin) / nucsPerPixel;

                drawArc(startLeftPix, startRightPix, endLeftPix, endRightPix, trackRectangle, context, feature.color);
                arcCount++;

                //System.out.println("    leftStart = " + leftStartNucPix);
                //System.out.println("    leftEnd = " + leftEndNucPix);
                //System.out.println("    rightStart = " + rightStartNucPix);
                //System.out.println("    rightEnd = " + rightEndNucPix);
                //System.out.println("    arcWidth = " + arcWidthPix);

                //drawArc(10, 210, 50, trackRectangle, context, ARC_COLOR_B);
                //drawArc(300, 500, 50, trackRectangle, context, ARC_COLOR_A);


                //System.out.println("Drew "+arcCount+" arcs");
            }
        }
        //System.out.println("");

        //draw a central horizontal line
        //Graphics2D g2D = context.getGraphic2DForColor(COLOR_CENTERLINE);

        //g2D.drawLine((int) trackRectangle.getX(), y,
        //        (int) trackRectangle.getMaxX(), y);
    }




    /**
     * Draw a filled arc between two regions of equal length
     *
     * @param startLeft  the starting position of the feature, whether on-screen or not
     * @param endLeft    the ending position of the feature, whether on-screen or not
     * @param pixelWidth  the width of the arc in pixels
     * @param trackRectangle
     * @param context
     * @param featureColor       the color specified for this feature.  May be null.
     */
    protected void drawArc(double startLeft, double startRight, double endLeft, double endRight,
                           Rectangle trackRectangle, RenderContext context, Color featureColor) {

        Color color;
        if (featureColor != null) {
            color = featureColor;
        } else {
            color = ARC_COLOR_A;
        }

        Graphics2D g2D = context.getGraphic2DForColor(color);
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
            g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        }
        //Height of top of an arc of maximum depth
        int maxPossibleArcHeight = (trackRectangle.height - 1) / 2;

        // Equation by G. Adam Stanislav from http://www.whizkidtech.redprince.net/bezier/circle/
        double handleLengthFactor = 4f*((double)Math.sqrt(2f)-1f)/3f;

        double outerRadius = (double) (endRight-startLeft)/2.0;
        double innerRadius = (double) (endLeft-startRight)/2.0;

        double arcWidth = Math.max(1.0, startRight-startLeft);

        int y = 0;
        if (dir>0){
            y = (int) trackRectangle.getMaxY();
        } else {
            y = (int) trackRectangle.getMinY();
        }

        // Define all control points
        // Use a minimum arc width of 1 pixel
        int outerLX = (int) (trackRectangle.getX() + startLeft);
        int outerLY = y;

        int outerLC1X = (int) (trackRectangle.getX() + outerLX);
        int outerLC1Y = (int) (outerLY - dir * handleLengthFactor*outerRadius);

        int outerLC2X = (int) (trackRectangle.getX() + outerLX+outerRadius - handleLengthFactor*outerRadius);
        int outerLC2Y = (int) (outerLY - dir * outerRadius);

        int outerCenterX = (int) (trackRectangle.getX() + outerLX+outerRadius);
        int outerCenterY = (int) (outerLY - dir * outerRadius);

        int outerRC1X = (int) (trackRectangle.getX() + outerLX+outerRadius + handleLengthFactor*outerRadius);
        int outerRC1Y = (int) (outerLY - dir * outerRadius);

        int outerRC2X = (int) (trackRectangle.getX() + endRight);
        int outerRC2Y = (int) (outerLY - dir * handleLengthFactor*outerRadius);

        int outerRX = (int) (trackRectangle.getX() + endRight);
        int outerRY = outerLY;

        int innerRX = (int) (trackRectangle.getX() + endRight - arcWidth);
        int innerRY = outerLY;

        int innerRC1X = (int) (trackRectangle.getX() + innerRX);
        int innerRC1Y = (int) (outerLY - dir * handleLengthFactor*innerRadius);

        int innerRC2X = (int) (trackRectangle.getX() + outerLX + outerRadius + handleLengthFactor*innerRadius);
        int innerRC2Y = (int) (outerCenterY + dir * arcWidth);

        int innerCenterX = (int) (trackRectangle.getX() + outerLX + outerRadius);
        int innerCenterY = (int) (outerCenterY + dir * arcWidth);

        int innerLC1X = (int) (trackRectangle.getX() + outerLX + outerRadius - handleLengthFactor*innerRadius);
        int innerLC1Y = (int) (outerCenterY + dir * arcWidth);

        int innerLC2X = (int) (trackRectangle.getX() + startLeft + arcWidth);
        int innerLC2Y = (int) (outerLY - dir * handleLengthFactor*innerRadius);

        int innerLX = (int) (trackRectangle.getX() + startLeft + arcWidth);
        int innerLY = outerLY;

        if (false){
            GeneralPath arcCage = new GeneralPath();
            arcCage.moveTo(outerLX, outerLY); // outer left
            arcCage.lineTo(outerLC1X, outerLC1Y);
            arcCage.lineTo(outerLC2X, outerLC2Y);
            arcCage.lineTo(outerCenterX, outerCenterY); // outer left control 1, outer left control 2, outer center
            arcCage.lineTo(outerRC1X, outerRC1Y);
            arcCage.lineTo(outerRC2X, outerRC2Y);
            arcCage.lineTo(outerRX, outerRY);// outer right control 1, outer right control 2, outer right
            arcCage.lineTo(innerRX, innerRY); // inner right
            arcCage.lineTo(innerRC1X, innerRC1Y);
            arcCage.lineTo(innerRC2X, innerRC2Y);
            arcCage.lineTo(innerCenterX, innerCenterY); // inner right control 1, inner right control 2, inner center
            arcCage.lineTo(innerLC1X, innerLC1Y);
            arcCage.lineTo(innerLC2X, innerLC2Y);
            arcCage.lineTo(innerLX, innerLY); // inner left control 1, inner left control 2, inner left
            arcCage.lineTo(outerLX, outerLY); // outer left
            arcCage.moveTo(outerLX, outerLY);
            arcCage.closePath();
            // Draw the shape
            //g2D.draw(arcCage);
            g2D.fill(arcCage);
        }

        // Create path
        if (true){
            GeneralPath arcPath = new GeneralPath();
            arcPath.moveTo(outerLX, outerLY); // outer left
            arcPath.curveTo(outerLC1X, outerLC1Y,
                    outerLC2X, outerLC2Y,
                    outerCenterX, outerCenterY); // outer left control 1, outer left control 2, outer center
            arcPath.curveTo(outerRC1X, outerRC1Y,
                    outerRC2X, outerRC2Y,
                    outerRX, outerRY);// outer right control 1, outer right control 2, outer right
            arcPath.lineTo(innerRX, innerRY); // inner right
            arcPath.curveTo(innerRC1X, innerRC1Y,
                    innerRC2X, innerRC2Y,
                    innerCenterX, innerCenterY); // inner right control 1, inner right control 2, inner center
            arcPath.curveTo(innerLC1X, innerLC1Y,
                    innerLC2X, innerLC2Y,
                    innerLX, innerLY); // inner left control 1, inner left control 2, inner left
            arcPath.lineTo(outerLX, outerLY); // outer left
            arcPath.moveTo(outerLX, outerLY);
            arcPath.closePath();

            // Draw the arc face
            g2D.fill(arcPath);
        }

        g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_DEFAULT);
        g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_DEFAULT);
    }
}