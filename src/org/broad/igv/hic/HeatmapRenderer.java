package org.broad.igv.hic;

import org.apache.commons.math.linear.RealMatrix;
import org.broad.igv.hic.data.Block;
import org.broad.igv.hic.data.ContactRecord;
import org.broad.igv.hic.data.DensityFunction;
import org.broad.igv.hic.data.MatrixZoomData;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.ui.Main;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Aug 11, 2010
 */
public class HeatmapRenderer {

    // TODO -- introduce a "model" in lieu of MainWindow pointer
    MainWindow mainWindow;

    public HeatmapRenderer(MainWindow mainWindow) {
        this.mainWindow = mainWindow;
    }

    public void render(int originX,
                       int originY,
                       int width,
                       int height,
                       MatrixZoomData zd,
                       Graphics g) {

        int chr1 = zd.getChr1();
        int chr2 = zd.getChr2();

        int maxX = originX + width;
        int maxY = originY + height;

        int x = originX;
        int y = originY;

        boolean isWholeGenome = zd.getChr1() == 0 && zd.getChr2() == 0;
        boolean sameChr = (chr1 == chr2);
        double binSizeMB = zd.getBinSize() / (isWholeGenome ? 1000.0 : 1000000.0);

        if (sameChr) {
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

        DensityFunction df = null;
        if (mainWindow.getDisplayOption() == MainWindow.DisplayOption.OE ||
                mainWindow.getDisplayOption() == MainWindow.DisplayOption.PEARSON) {
            df = mainWindow.getDensityFunction(zd.getZoom());
        }


        // Iterate through blocks overlapping visible region
        ColorScale colorScale = mainWindow.getColorScale();

        MainWindow.DisplayOption displayOption = mainWindow.getDisplayOption();
        if (displayOption == MainWindow.DisplayOption.PEARSON && df != null) {
            RealMatrix rm = zd.getPearsons(df);
            renderMatrix(originX, originY, zd.getBinSize(), rm, colorScale, g);

        } else {
            List<Block> blocks = zd.getBlocksOverlapping(x, y, maxX, maxY);
            for (Block b : blocks) {
                renderBlock(originX, originY, chr1, chr2, binSizeMB, b, colorScale, df, g);
            }
        }
    }

    private void renderBlock(int originX, int originY, int chr1, int chr2, double binSizeMB, Block b,
                             ColorScale colorScale, DensityFunction df, Graphics g) {

        MainWindow.DisplayOption displayOption = mainWindow.getDisplayOption();
        double binSizeMB2 = binSizeMB * binSizeMB;
        boolean sameChr = (chr1 == chr2);

        ContactRecord[] recs = b.getContactRecords();
        if (recs != null) {
            for (int i = 0; i < recs.length; i++) {
                ContactRecord rec = recs[i];

                Color color = null;
                double score;
                if (displayOption == MainWindow.DisplayOption.OE && df != null) {
                    int x = rec.getX();// * binSize;
                    int y = rec.getY();// * binSize;
                    int dist = Math.abs(x - y);
                    double expected = df.getDensity(chr1, dist);
                    score = rec.getCounts() / expected;
                    score = Math.log10(score);
                } else {
                    score = rec.getCounts() / binSizeMB2;
                }

                color = colorScale.getColor((float) score);
                int px = (rec.getX() - originX);
                int py = (rec.getY() - originY);

                g.setColor(color);
                g.fillRect(px, py, MainWindow.BIN_PIXEL_WIDTH, MainWindow.BIN_PIXEL_WIDTH);

                if (sameChr && (rec.getX() != rec.getY())) {
                    px = (rec.getY() - originX);
                    py = (rec.getX() - originY);
                    g.fillRect(px, py, MainWindow.BIN_PIXEL_WIDTH, MainWindow.BIN_PIXEL_WIDTH);
                }
            }
        }
    }

    private void renderMatrix(int originX, int originY, int binSize, RealMatrix rm,
                              ColorScale colorScale, Graphics g) {

        MainWindow.DisplayOption displayOption = mainWindow.getDisplayOption();

        int nBinsX = rm.getColumnDimension();
        int nBinsY = rm.getRowDimension();

        for (int i = 0; i < nBinsX; i++) {
            for (int j = 0; j < nBinsY; j++) {
                double score = rm.getEntry(i, j);
                //float logScore = (float) Math.log10(score);
                Color color = score == 0 ? Color.black : colorScale.getColor((float) score);
                int px = i - originX;
                int py = j - originY;
                g.setColor(color);
                g.fillRect(px, py, MainWindow.BIN_PIXEL_WIDTH, MainWindow.BIN_PIXEL_WIDTH);
                // Assuming same chromosome
                if (i != j) {
                    px = (j - originX);
                    py = (i - originY);
                    g.fillRect(px, py, MainWindow.BIN_PIXEL_WIDTH, MainWindow.BIN_PIXEL_WIDTH);
                }
            }
        }
    }


}
