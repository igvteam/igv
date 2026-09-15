package org.igv.alignment.fiberseq;

import org.igv.alignment.Alignment;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;

import java.awt.*;
import java.util.List;

import static org.igv.prefs.Constants.*;

/**
 * Draws fibertools nucleosomes and methylation sensitive patches (MSPs) over an alignment.  MSPs with a non-zero FIRE
 * quality are drawn in the FIRE color, more opaque as the quality increases.  Colors are user preferences.
 */
public class FiberseqRenderer {

    private static Color fireColor;
    private static Color[] fireColors;

    public static void draw(Alignment alignment, double bpStart, double locScale, Rectangle rowRect, Graphics g,
                            boolean leaveMargin) {
        FiberseqAnnotations annotations = alignment.getFiberseqAnnotations();
        if (annotations == null) {
            return;
        }
        int h = Math.max(1, rowRect.height - (leaveMargin ? 2 : 0));
        int y = rowRect.y;

        IGVPreferences prefs = PreferencesManager.getPreferences();
        g.setColor(prefs.getAsColor(FIBERSEQ_NUCLEOSOME_COLOR));
        drawIntervals(annotations.getNucleosomes(), null, null, bpStart, locScale, rowRect, g, y, h);
        drawIntervals(annotations.getMsps(), prefs.getAsColor(FIBERSEQ_MSP_COLOR),
                getFireColors(prefs.getAsColor(FIBERSEQ_FIRE_COLOR)), bpStart, locScale, rowRect, g, y, h);
    }

    private static synchronized Color[] getFireColors(Color color) {
        if (!color.equals(fireColor)) {
            Color[] colors = new Color[256];
            for (int q = 0; q < colors.length; q++) {
                colors[q] = new Color(color.getRed(), color.getGreen(), color.getBlue(), 100 + (155 * q) / 255);
            }
            fireColors = colors;
            fireColor = color;
        }
        return fireColors;
    }

    /**
     * Draw intervals in the current color, or for MSPs (mspColor non-null) in the MSP or FIRE color.
     */
    private static void drawIntervals(List<FiberseqAnnotations.Interval> intervals, Color mspColor,
                                      Color[] fireColors, double bpStart, double locScale, Rectangle rowRect,
                                      Graphics g, int y, int h) {
        for (FiberseqAnnotations.Interval interval : intervals) {
            int pStart = (int) ((interval.start() - bpStart) / locScale);
            int pEnd = (int) ((interval.end() - bpStart) / locScale);
            // Intervals are in molecular order, which is descending reference order for reverse-strand reads,
            // so an interval past either edge of the view does not end the loop
            if (pEnd < rowRect.x || pStart > rowRect.getMaxX()) {
                continue;
            }
            if (mspColor != null) {
                g.setColor(interval.quality() > 0 ? fireColors[Math.min(255, interval.quality())] : mspColor);
            }
            g.fillRect(pStart, y, Math.max(1, pEnd - pStart), h);
        }
    }
}
