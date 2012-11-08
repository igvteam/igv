package org.broad.igv.hic.track;

import org.broad.igv.hic.Context;
import org.broad.igv.hic.HiC;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.track.RenderContext;

import java.awt.*;

/**
 * @author jrobinso
 *         Date: 11/2/12
 *         Time: 9:41 AM
 */
public class HiCDataTrack extends HiCTrack {

    static int TRACK_MARGIN = 2;
    HiC hic;
    HiCDataAdapter dataSource;
    boolean logScale = false;


    public HiCDataTrack(HiC hic, HiCDataAdapter dataSource) {
        this.hic = hic;
        this.dataSource = dataSource;
        this.logScale = dataSource.isLogScale();

    }

    @Override
    public void render(Graphics2D g2d, Context context, Rectangle rect) {


        int h = rect.height;

        String chr = context.getChromosome().getName();
        int start = context.getBinOrigin();
        int end = start + (int) (rect.width / context.getScaleFactor());
        HiCDataAdapter.WeightedSum[] data = dataSource.getData(chr, start, end);

        Color posColor = dataSource.getColor();
        Color negColor = dataSource.getAltColor();

        // Get the Y axis definition, consisting of minimum, maximum, and base value.  Often
        // the base value is == min value which is == 0.

        DataRange dataRange = dataSource.getDataRange();
        float maxValue = dataRange.getMaximum();
        float baseValue = dataRange.getBaseline();
        float minValue = dataRange.getMinimum();
        boolean isLog = dataRange.isLog();

        if (isLog) {
            minValue = (float) (minValue == 0 ? 0 : Math.log10(minValue));
            maxValue = (float) Math.log10(maxValue);
        }


        // Calculate the Y scale factor.

        double delta = (maxValue - minValue);
        double yScaleFactor = (rect.height - TRACK_MARGIN) / delta;

        // Calculate the Y position in pixels of the base value.  Clip to bounds of rectangle
        double baseDelta = maxValue - baseValue;
        int baseY = (int) (rect.y + baseDelta * yScaleFactor);
        if (baseY < rect.y) {
            baseY = rect.y;
        } else if (baseY > rect.y + (rect.height - TRACK_MARGIN)) {
            baseY = rect.y + (rect.height - TRACK_MARGIN);
        }

        for (int i = 0; i < data.length; i++) {

            HiCDataAdapter.WeightedSum d = data[i];
            if (d == null) continue;

            int bin = d.getBinNumber() - start;
            int xPixelLeft = rect.x + (int) ((bin) * context.getScaleFactor()); //context.getScreenPosition (genomicPosition);
            int xPixelRight = rect.x + (int) ((bin + 1) * context.getScaleFactor());
            int dx = xPixelRight - xPixelLeft;

            if (xPixelRight < rect.x) {
                continue;
            } else if (xPixelLeft > rect.x + rect.width) {
                break;
            }

            double dataY = d.getValue();
            if (isLog && dataY <= 0) {
                continue;
            }

            if (!Double.isNaN(dataY)) {

                // Compute the pixel y location.  Clip to bounds of rectangle.
                double dy = isLog ? Math.log10(dataY) - baseValue : (dataY - baseValue);
                int pY = baseY - (int) (dy * yScaleFactor);
                if (pY < rect.y) {
                    pY = rect.y;
                } else if (pY > rect.y + (rect.height - TRACK_MARGIN)) {
                    pY = rect.y + (rect.height - TRACK_MARGIN);
                }

                Color color = (dataY >= baseValue) ? posColor : negColor;
                g2d.setColor(color);

                if (dx <= 1) {
                    g2d.drawLine(xPixelRight, baseY, xPixelRight, pY);
                } else {
                    if (pY > baseY) {
                        g2d.fillRect(xPixelRight, baseY, dx, pY - baseY);

                    } else {
                        g2d.fillRect(xPixelRight, pY, dx, baseY - pY);
                    }
                }
            }
        }
    }

    public String getName() {
        return dataSource.getName();
    }

    @Override
    public Color getColor() {
        return dataSource.getColor();
    }
}
