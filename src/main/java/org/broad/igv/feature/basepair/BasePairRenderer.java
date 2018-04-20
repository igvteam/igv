package org.broad.igv.feature.basepair;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.basepair.BasePairTrack.*;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.awt.geom.GeneralPath;

// TODO: add alternative visualization for very zoomed-out views

public class BasePairRenderer {

    private static Logger log = Logger.getLogger(BasePairRenderer.class);

    public void draw(BasePairData data,
                     RenderContext context,
                     Rectangle trackRectangle,
                     RenderOptions renderOptions){

        double nucsPerPixel = context.getScale();
        double origin = context.getOrigin();
        //The location of the first base that is loaded
        int start = Math.max(0, (int) origin - 1);
        //The location of the last base that is loaded
        int end = (int) (origin + trackRectangle.width * nucsPerPixel) + 1;
        if (end <= start) return;

        java.util.List<BasePairFeature> featureList = data.getFeatures(context.getChr());

        if (featureList != null) {

            // render in order of indexed colors
            // FIXME: this is inefficient/slow (but not too bad in practice for most current use cases)
            for (int colorIndex=0; colorIndex<renderOptions.getColors().size(); ++colorIndex) {
                for (BasePairFeature feature : featureList) {

                    if(feature.startLeft > context.getEndLocation()) break;
                    else if(feature.endRight < context.getOrigin()) continue;
                    else if (feature.colorIndex != colorIndex) continue;

                    double startLeftPix = (feature.startLeft - origin) / nucsPerPixel;
                    double startRightPix = (feature.startRight + 1.0 - origin) / nucsPerPixel;
                    double endLeftPix = (feature.endLeft - origin) / nucsPerPixel;
                    double endRightPix = (feature.endRight + 1.0 - origin) / nucsPerPixel;

                    Color color = ColorUtilities.stringToColor(renderOptions.getColors().get(feature.colorIndex));

                    boolean drawOutline=false;
                    drawArc(startLeftPix, startRightPix, endLeftPix, endRightPix,
                            trackRectangle, context, color,
                            renderOptions.getArcDirection(), drawOutline);

                }
            }
        }
    }




    /**
     * Draw a filled arc between two regions of equal length
     *
     * @param startLeft  the starting position of the feature, whether on-screen or not (outer left corner of arc)
     * @param startRight (inner left corner of arc)
     * @param endLeft    the ending position of the feature, whether on-screen or not (inner right corner of arc)
     * @param endRight   (outer right corner of arc)
     * @param trackRectangle
     * @param context
     * @param featureColor       the color specified for this feature.  May be null.
     * @param arcDirection up or down
     * @param drawOutline render outline around arc shape so always visible even when zoomed out
     */
    protected void drawArc(double startLeft, double startRight, double endLeft, double endRight,
                           Rectangle trackRectangle, RenderContext context, Color featureColor,
                           ArcDirection arcDirection, boolean drawOutline) {

        Color color;
        if (featureColor != null) {
            color = featureColor;
        } else {
            color = new Color(50, 50, 150, 140);
        }

        Graphics2D g2D = context.getGraphic2DForColor(color);

        //Height of top of an arc of maximum depth
        int maxPossibleArcHeight = (trackRectangle.height - 1) / 2;

        // Equation by G. Adam Stanislav from http://www.whizkidtech.redprince.net/bezier/circle/
        double handleLengthFactor = 4f*((double)Math.sqrt(2f)-1f)/3f;

        double outerRadius = (double) (endRight-startLeft)/2.0;
        double innerRadius = (double) (endLeft-startRight)/2.0;

        double arcWidth = Math.max(1.0, startRight-startLeft);

        int dir  = 1; // 1 for up, -1 for down
        if (arcDirection == ArcDirection.DOWN) dir = -1;

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

        // Draw outline so thin arcs don't vanish when zoomed out
        if (drawOutline) g2D.draw(arcPath);
        // Draw the arc face
        g2D.fill(arcPath);

        g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_DEFAULT);
        g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_DEFAULT);
    }
}
