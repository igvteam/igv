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

    public void render(int originX, int originY, MatrixZoomData zd, int binWidth, double maxCount, Graphics g, Rectangle bounds, Color background) {

        int chr1 = zd.getChr1(); //Context.chr1Index;
        int chr2 = zd.getChr2(); //Context.chr2Index;


        final double boundsRight = bounds.getMaxX();
        final double boundsLeft = bounds.getX();
        final double boundsLower = bounds.getMaxY();
        final double boundsUpper = bounds.getY();

        // TODO Lookup zd object from chr1Index, chr2Index, zoom

        // Iterate through blocks overlapping visible region
        int binSize = zd.getBinSize();
        int maxX = originX + (int) (bounds.getWidth() * binSize / binWidth) + 1;
        int maxY = originY + (int) (bounds.getHeight() * binSize / binWidth) + 1;
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

        for (Block b : blocks) {

            ContactRecord[] recs = b.getContactRecords();
            if (recs != null) {
                for (int i = 0; i < recs.length; i++) {
                    ContactRecord rec = recs[i];

                    Color color = null;
                    double binSizeHKB = binSize / 100000.0;
                    double score = rec.getCounts() / (binSizeHKB * binSizeHKB);
                    if (maxCount > 0 && score > 2 * maxCount) {
                        color = Color.ORANGE;
                    } else {
                        float alpha = (float) Math.max(0.05f, Math.min(1.0f, score / maxCount));
                        color = getColor(alpha, background);
                    }

                    int px = bounds.x + (rec.getX() - originX / binSize) * binWidth;
                    int py = bounds.y + (rec.getY() - originY / binSize) * binWidth;


                    g.setColor(color);
                    if (px <= boundsRight && (px + binWidth) >= boundsLeft &&
                            py <= boundsLower && py + binWidth >= boundsUpper) {
                        g.fillRect(px, py, binWidth, binWidth);
                    }


                    if (chr1 == chr2) {
                        px = bounds.x + (rec.getY() - originX / binSize) * binWidth;
                        py = bounds.y + (rec.getX() - originY / binSize) * binWidth;
                        if (px <= boundsRight && (px + binWidth) >= boundsLeft &&
                                py <= boundsLower && py + binWidth >= boundsUpper) {
                            g.fillRect(px, py, binWidth, binWidth);
                        }

                    }
                }

            }
        }
    }

    // Cache colors,  100 distinct colors should be enough

    static Map<Integer, Color> colorCache = new Hashtable();

    static Color getColor(float alpha, Color background) {

        float [] comps = background.getColorComponents(null);

        int idx = (int) (100 * alpha);
        Color c = colorCache.get(idx);
        if (c == null) {
            float rAlpha = Math.max(0.05f, Math.min(1.0f, 0.01f * idx));
            float red =  ((1-rAlpha) * comps[0] + rAlpha);
            float green =  ((1-rAlpha) * comps[1]);
            float blue = ((1-rAlpha) * comps[2]);
            c = new Color(red, green, blue);
            colorCache.put(idx, c);
        }
        return c;
    }

}
