package org.broad.igv.ui.supdiagram;

import org.broad.igv.feature.Strand;

import java.awt.*;
import java.awt.geom.Point2D;

/**
 * An arrow shape that represents  a single Supplementary Read in the supplementary alignment diagram
 */
public class AlignmentArrow extends Polygon {
    /** The width of the arrow tip in pixels */
    public static final int ARROW_PX_WIDTH = 5;
    final Strand strand;

    public AlignmentArrow(int midline, int height, int left, int right, Strand strand) {
        super();
        final int floor = (midline - height);
        final int startAdjusted = left - (strand == Strand.NEGATIVE ? ARROW_PX_WIDTH : 0);
        final int endAdjusted = right + (strand == Strand.POSITIVE ? ARROW_PX_WIDTH : 0);
         /*
                1       2
             0 <|=======|> 3
                5       4
        */
        final int[] xPoly = {startAdjusted, left, right, endAdjusted, right, left};
        final int[] yPoly = {floor + height / 2, floor, floor, floor + height / 2, floor + height, floor + height};
        this.xpoints = xPoly;
        this.ypoints = yPoly;
        this.npoints = xPoly.length;
        this.strand = strand;
        invalidate();
    }

    public Point2D getTip() {
        int x = strand == Strand.NEGATIVE ? xpoints[0] : xpoints[3];
        return new Point2D.Double(x, ypoints[0]);
    }

    public Point2D getTail() {
        int x = strand == Strand.NEGATIVE ? xpoints[3] : xpoints[0];
        return new Point2D.Double(x, ypoints[0]);
    }
}
