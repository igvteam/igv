package org.igv.feature.mut;

import org.igv.feature.LocusScore;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.Renderer;
import org.igv.track.DataType;
import org.igv.track.RenderContext;
import org.igv.track.Track;

import java.awt.*;
import java.util.List;

public class MutRenderer implements Renderer<Mutation> {

    @Override
    public void render(List<Mutation> features, RenderContext context, Rectangle rect, Track track) {

        //AbstractColorScale colorScale = track.getColorScale();
        double origin = context.getOrigin();
        double locScale = context.getScale();

        double maxX = rect.getMaxX();
        int minY = (int) rect.getMinY();
        int height = (int) rect.getHeight();
        int lastPEnd = 0;
        int lastPStart = 0;
        int lastW = 0;

        for (LocusScore score : features) {
            if (lastPStart > maxX) {
                break;
            }

            // Note -- don't cast these to an int until the range is checked,
            // otherwise could get an overflow.
            float fStart = (float) ((score.getStart() - origin) / locScale);
            float fEnd = (float) ((score.getEnd() - origin) / locScale);
            // float fw = fEnd - fStart;
            int pStart = (int) fStart;
            int pEnd = (int) fEnd;

            int min = 1;
            int w = Math.max(min, pEnd - pStart);

            Color graphColor = Color.blue; //colorScale.getColor(track.logScaleData(score.getScore()));

            if ((pStart + w) >= 0 && (lastPStart <= maxX)) {

                // Minimum width for DNA methylation
                if (track.getDataType() == DataType.DNA_METHYLATION) {
                    if (w < 6) {
                        int pMid = (pStart + pEnd) / 2;
                        pStart = pMid - 3;
                        w = 6;
                    }
                }

                // This test handles the rather pathological case where the previous feature was 1 pixel wide, and
                // the current feature overlaps it because of roundoff error when scaling.
                else if (pStart < lastPEnd && w > 1 && lastW == 1) {
                    pStart++;
                    w--;
                }

                Graphics2D g2D = context.getGraphic2DForColor(graphColor);
                if (pStart < maxX) {
                    // Clip at edges
                    int pLeft = Math.max(rect.x, pStart);
                    int pRight = Math.min(rect.x + rect.width, pStart + w);
                    int adjustedW = pRight - pLeft;
                    g2D.fillRect(pLeft, minY, adjustedW, height);
                }
            }
            lastPStart = pStart;
            lastPEnd = pStart + w;
            lastW = w;
        }
    }
}
