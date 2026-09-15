package org.igv.alignment.fiberseq;

import org.igv.alignment.Alignment;

import java.awt.*;
import java.util.List;

/**
 * Draws fibertools nucleosomes (gray) and methylation sensitive patches (purple) over an alignment.  MSPs with a
 * non-zero FIRE quality are drawn red, more opaque as the quality increases.
 */
public class FiberseqRenderer {

    static final Color NUCLEOSOME_COLOR = new Color(150, 150, 150);
    static final Color MSP_COLOR = new Color(147, 112, 219);
    private static final Color[] FIRE_COLORS = new Color[256];

    static {
        for (int q = 0; q < FIRE_COLORS.length; q++) {
            FIRE_COLORS[q] = new Color(200, 0, 0, 100 + (155 * q) / 255);
        }
    }

    public static void draw(Alignment alignment, double bpStart, double locScale, Rectangle rowRect, Graphics g,
                            boolean leaveMargin) {
        FiberseqAnnotations annotations = alignment.getFiberseqAnnotations();
        if (annotations == null) {
            return;
        }
        int h = Math.max(1, rowRect.height - (leaveMargin ? 2 : 0));
        int y = rowRect.y;

        g.setColor(NUCLEOSOME_COLOR);
        drawIntervals(annotations.getNucleosomes(), false, bpStart, locScale, rowRect, g, y, h);
        drawIntervals(annotations.getMsps(), true, bpStart, locScale, rowRect, g, y, h);
    }

    private static void drawIntervals(List<FiberseqAnnotations.Interval> intervals, boolean msp, double bpStart,
                                      double locScale, Rectangle rowRect, Graphics g, int y, int h) {
        for (FiberseqAnnotations.Interval interval : intervals) {
            int pStart = (int) ((interval.start() - bpStart) / locScale);
            int pEnd = (int) ((interval.end() - bpStart) / locScale);
            // Intervals are in molecular order, which is descending reference order for reverse-strand reads,
            // so an interval past either edge of the view does not end the loop
            if (pEnd < rowRect.x || pStart > rowRect.getMaxX()) {
                continue;
            }
            if (msp) {
                g.setColor(interval.quality() > 0 ? FIRE_COLORS[Math.min(255, interval.quality())] : MSP_COLOR);
            }
            g.fillRect(pStart, y, Math.max(1, pEnd - pStart), h);
        }
    }
}
