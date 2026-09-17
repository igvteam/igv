package org.igv.alignment.fiberseq;

import org.igv.alignment.Alignment;
import org.igv.alignment.Gap;
import org.igv.alignment.SAMAlignment;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;

import java.awt.*;
import java.util.List;

import static org.igv.prefs.Constants.*;

/**
 * Draws fibertools nucleosomes and methylation sensitive patches (MSPs) over an alignment.  MSPs with a non-zero FIRE
 * quality are drawn in the FIRE color, more opaque as the quality increases.  Colors are user preferences.
 * <p>
 * Where an interval spans a gap (deletion or skipped region) it is drawn as a line matching the gap line, so the gap
 * stays visible.
 */
public class FiberseqRenderer {

    private static Color fireColor;
    private static Color[] fireColors;

    /**
     * @param minGapWidth gaps narrower than this are not drawn as gaps (small indels hidden), so the overlay is solid
     *                    across them
     */
    public static void draw(Alignment alignment, double bpStart, double locScale, Rectangle rowRect, Graphics g,
                            boolean leaveMargin, int minGapWidth) {
        MolecularAnnotations annotations = alignment.getMolecularAnnotations();
        if (annotations == null) {
            return;
        }
        int h = Math.max(1, rowRect.height - (leaveMargin ? 2 : 0));
        Geometry geom = new Geometry(bpStart, locScale, rowRect, g, rowRect.y, h, alignment.getGaps(), minGapWidth);

        IGVPreferences prefs = PreferencesManager.getPreferences();
        g.setColor(prefs.getAsColor(MA_NUCLEOSOME_COLOR));
        drawIntervals(annotations.getNucleosomes(), null, null, geom);
        drawIntervals(annotations.getMsps(), prefs.getAsColor(MA_MSP_COLOR),
                getFireColors(prefs.getAsColor(MA_FIRE_COLOR)), geom);
    }

    private record Geometry(double bpStart, double locScale, Rectangle rowRect, Graphics g, int y, int h,
                            List<Gap> gaps, int minGapWidth) {
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
    private static void drawIntervals(List<MolecularAnnotations.Interval> intervals, Color mspColor,
                                      Color[] fireColors, Geometry geom) {
        for (MolecularAnnotations.Interval interval : intervals) {
            int pStart = (int) ((interval.start() - geom.bpStart) / geom.locScale);
            int pEnd = (int) ((interval.end() - geom.bpStart) / geom.locScale);
            // Intervals are in molecular order, which is descending reference order for reverse-strand reads,
            // so an interval past either edge of the view does not end the loop
            if (pEnd < geom.rowRect.x || pStart > geom.rowRect.getMaxX()) {
                continue;
            }
            if (mspColor != null) {
                geom.g.setColor(interval.quality() > 0 ? fireColors[Math.min(255, interval.quality())] : mspColor);
            }
            drawInterval(interval.start(), interval.end(), geom);
        }
    }

    private static void drawInterval(int start, int end, Geometry geom) {
        int pos = start;
        if (geom.gaps != null) {
            for (Gap gap : geom.gaps) {
                int gapStart = gap.getStart();
                int gapEnd = gapStart + gap.getnBases();
                if (gapEnd <= pos || gap.getnBases() < geom.minGapWidth) {
                    continue;
                } else if (gapStart >= end) {
                    break;
                }
                fillSpan(pos, gapStart, geom);
                int lineEnd = Math.min(gapEnd, end);
                // Same geometry as the gap line: 2px through the row center for deletions in rows taller than 5px
                int thickness = (gap.getType() == SAMAlignment.DELETION && geom.h > 5) ? 2 : 1;
                int pLineStart = (int) ((Math.max(pos, gapStart) - geom.bpStart) / geom.locScale);
                int pLineEnd = (int) ((lineEnd - geom.bpStart) / geom.locScale);
                geom.g.fillRect(pLineStart, geom.y + geom.h / 2 - thickness / 2, Math.max(1, pLineEnd - pLineStart),
                        thickness);
                pos = lineEnd;
            }
        }
        fillSpan(pos, end, geom);
    }

    private static void fillSpan(int start, int end, Geometry geom) {
        if (start >= end) {
            return;
        }
        int pStart = (int) ((start - geom.bpStart) / geom.locScale);
        int pEnd = (int) ((end - geom.bpStart) / geom.locScale);
        geom.g.fillRect(pStart, geom.y, Math.max(1, pEnd - pStart), geom.h);
    }
}
