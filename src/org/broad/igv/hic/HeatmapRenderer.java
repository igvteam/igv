package org.broad.igv.hic;

import org.broad.igv.hic.data.Block;
import org.broad.igv.hic.data.ContactRecord;
import org.broad.igv.hic.data.MatrixZoomData;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Aug 11, 2010
 */
public class HeatmapRenderer {

    public void render(int originX,
                       int originY,
                       int width,
                       int height,
                       MatrixZoomData zd,
                       double maxCount,
                       Graphics g,
                       Color background) {

        int chr1 = zd.getChr1();
        int chr2 = zd.getChr2();


        int maxX = originX + width;
        int maxY = originY + height;

        // TODO Lookup zd object from chr1Index, chr2Index, zoom

        // Iterate through blocks overlapping visible region
        //int binSize = zd.getBinSize();
        int x = originX;
        int y = originY;

        if (chr1 == chr2) {
            // Data is transposable, transpose if neccessary.  Convention is to use lower diagonal
            if (x > y) {
                x = originY;
                y = originX;
            }
            if (maxX > maxY) {
                int tmp = maxX;
                maxX = maxY;
                maxY = tmp;
            }
        }


        List<Block> blocks = zd.getBlocksOverlapping(x, y, maxX, maxY);

        boolean isWholeGenome = zd.getChr1() == 0 && zd.getChr2() == 0;

        for (Block b : blocks) {

            ContactRecord[] recs = b.getContactRecords();
            if (recs != null) {
                for (int i = 0; i < recs.length; i++) {
                    ContactRecord rec = recs[i];

                    Color color = null;
                    double binSizeMB = zd.getBinSize() / (isWholeGenome ? 1000.0 :   1000000.0);
                    double score = rec.getCounts() / (binSizeMB * binSizeMB);
                    //if (maxCount > 0 && score > 2 * maxCount) {
                    //    color = Color.ORANGE;
                    //} else {
                        float alpha = (float) Math.max(0.05f, Math.min(1.0f, score / maxCount));
                        color = getColor(alpha, background);
                    //}

                    int px = (rec.getX() - originX );
                    int py = (rec.getY() - originY );


                    g.setColor(color);
                    g.fillRect(px, py, 1, 1);


                    if (chr1 == chr2) {
                        px = (rec.getY() - originX );
                        py =(rec.getX() - originY );
                        g.fillRect(px, py, 1, 1);
                    }


                }

            }
        }
    }

    // Cache colors,  100 distinct colors should be enough

    static Map<Integer, Color> colorCache = new Hashtable();

    static Color getColor(float alpha, Color background) {

        float[] comps = background.getColorComponents(null);

        int idx = (int) (100 * alpha);
        Color c = colorCache.get(idx);
        if (c == null) {
            float rAlpha = Math.max(0.05f, Math.min(1.0f, 0.01f * idx));
            float red = ((1 - rAlpha) * comps[0] + rAlpha);
            float green = ((1 - rAlpha) * comps[1]);
            float blue = ((1 - rAlpha) * comps[2]);
            c = new Color(red, green, blue);
            colorCache.put(idx, c);
        }
        return c;
    }

}
