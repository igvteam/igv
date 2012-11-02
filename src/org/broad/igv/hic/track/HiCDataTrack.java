package org.broad.igv.hic.track;

import org.broad.igv.hic.Context;
import org.broad.igv.hic.HiC;

import java.awt.*;

/**
 * @author jrobinso
 *         Date: 11/2/12
 *         Time: 9:41 AM
 */
public class HiCDataTrack extends HiCTrack {

    HiC hic;
    HiCDataAdapter dataSource;
    double dataMax;


    public HiCDataTrack(HiC hic, HiCDataAdapter dataSource) {
        this.hic= hic;
        this.dataSource = dataSource;
        this.dataMax = dataSource.getMax();

    }

    @Override
    public void render(Graphics2D g2d, Context context, Rectangle rect) {


        int h = rect.height;

       String chr = context.getChromosome().getName();
        int start = context.getBinOrigin();
        int end = start + (int) (rect.width / context.getScaleFactor());
        HiCDataAdapter.WeightedSum[] data = dataSource.getData(chr, start, end);


        for (int bin = 0; bin < data.length; bin++) {

            HiCDataAdapter.WeightedSum d = data[bin];
            if (d == null) continue;

            int xPixelLeft = rect.x + (int) ((bin) * context.getScaleFactor()); //context.getScreenPosition (genomicPosition);
            int xPixelRight = rect.x + (int) ((bin + 1) * context.getScaleFactor());


            if (xPixelRight < rect.x) {
                continue;
            } else if (xPixelLeft > rect.x + rect.width) {
                break;
            }

            System.out.println(dataMax);

            double mid = 0;
            double x = d.getValue() - mid;
            double max = dataMax - mid;

            int myh = (int) ((x / max) * h);
            //if (x > 0) {
            g2d.fillRect(xPixelLeft, rect.y + h - myh, (xPixelRight - xPixelLeft + 1), myh);
            //} else {
            //    g2d.fillRect(xPixelLeft, rect.y + h, (xPixelRight - xPixelLeft), -myh);
            //}

        }
    }
}
