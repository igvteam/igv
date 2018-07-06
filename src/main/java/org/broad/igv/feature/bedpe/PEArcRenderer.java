package org.broad.igv.feature.bedpe;

import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.awt.geom.GeneralPath;
import java.util.List;

import static org.broad.igv.feature.bedpe.PEArcRenderer.Direction.DOWN;
import static org.broad.igv.feature.bedpe.PEArcRenderer.Direction.UP;


/**
 * Created by jrobinso on 6/29/18.
 */
public class PEArcRenderer {

    enum Direction {UP, DOWN}

    double theta = Math.toRadians(45);
    double sinTheta = Math.sin(theta);
    double cosTheta = Math.cos(theta);
    Direction direction = DOWN;
    int lineThickness = 1;


    public void render(List<BedPEFeature> featureList, RenderContext context, Rectangle trackRectangle, Track track) {

        double origin = context.getOrigin();
        double locScale = context.getScale();


        for (BedPEFeature feature : featureList) {

            // Note -- don't cast these to an int until the range is checked.
            // could get an overflow.
            double pixelStart = ((feature.getStart() - origin) / locScale);
            double pixelEnd = ((feature.getEnd() - origin) / locScale);

            // If the any part of the feature fits in the Track rectangle draw it
            if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {

                int s = (feature.start1 + feature.end1) / 2;
                int e = (feature.start2 + feature.end2) / 2;

                double ps = (s - origin) / locScale;
                double pe = (e - origin) / locScale;
                drawArc(context, trackRectangle, track, ps, pe);

            }
        }
    }

    private void drawArc(RenderContext context, Rectangle trackRectangle, Track track, double x1, double x2) {

        double pixelStart = Math.min(x1, x2);
        double pixelEnd = Math.max(x1, x2);

        Color color = track.getColor();
        Graphics2D g = context.getGraphic2DForColor(color);
        if(lineThickness > 1) g.setStroke(new BasicStroke(lineThickness));

        int w = (int) (pixelEnd - pixelStart);
        if (w < 3) {
            w = 3;
            pixelStart--;
        }

        double a = w / 2;
        double r = a / sinTheta;
        double b = cosTheta * r;
        double c = r - b;
        double d = r - a;
        double x = pixelStart - d;
        final int trackBaseLine = trackRectangle.y + trackRectangle.height;
        double y = trackBaseLine - c;

        double angleSt = 90 - Math.toDegrees(theta);
        double ext = Math.toDegrees(2 * theta);

        Arc2D.Double arcPath;
        if (direction == UP) {
            arcPath = new Arc2D.Double(x, y, 2 * r, 2 * r, angleSt, ext, Arc2D.OPEN);
        } else {
            y = trackRectangle.y - 2 * r + c;
            angleSt = 180 + angleSt;
            arcPath = new Arc2D.Double(x, y, 2 * r, 2 * r, angleSt, ext, Arc2D.OPEN);
        }

        g.draw(arcPath);
    }




    static void sagittaCalcs() {

        for(double theta = 0; theta <= Math.PI / 2;  theta += Math.PI/100) {

            double num = 1 - Math.cos(theta);
            double denom = Math.sin(theta);

            System.out.println(theta + "\t" + theta*Math.PI/2 + "\t" + (num / denom));
        }


    }

}
