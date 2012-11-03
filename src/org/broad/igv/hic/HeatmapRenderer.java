package org.broad.igv.hic;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.hic.data.Block;
import org.broad.igv.hic.data.ContactRecord;
import org.broad.igv.hic.data.DensityFunction;
import org.broad.igv.hic.data.MatrixZoomData;
import org.broad.igv.hic.matrix.BasicMatrix;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.util.collections.DoubleArrayList;

import java.awt.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 * @date Aug 11, 2010
 */
public class HeatmapRenderer {

    // TODO -- introduce a "model" in lieu of MainWindow pointer
    HiC hic;
    MainWindow mainWindow;

    private ContinuousColorScale observedColorScale;
    private ColorScale oeColorScale;
    private ColorScale pearsonColorScale;

    Map<MatrixZoomData, ContinuousColorScale> observedColorScaleMap = new HashMap<MatrixZoomData, ContinuousColorScale>();


    public HeatmapRenderer(MainWindow mainWindow, HiC hic) {
        this.mainWindow = mainWindow;
        this.hic = hic;

        oeColorScale = new HiCColorScale();
        pearsonColorScale = new HiCColorScale();
    }

    public void render(int originX,
                       int originY,
                       int width,
                       int height,
                       final MatrixZoomData zd,
                       MainWindow.DisplayOption displayOption,
                       Graphics2D g) {

        g.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_SPEED);
        g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);

        int chr1 = zd.getChr1();
        int chr2 = zd.getChr2();


        int x = originX;
        int y = originY;

        boolean isWholeGenome = zd.getChr1() == 0 && zd.getChr2() == 0;
        boolean sameChr = (chr1 == chr2);
        int zoomMultiplier = zd.getZoomMultiplier();
        double mb;
        if (isWholeGenome) {
            mb = (double) zoomMultiplier * 5;
        } else {
            mb = (double) (zoomMultiplier) / 200;
        }
        double colorScaleFactor = mb * mb;


        if (sameChr) {
            // Data is transposable, transpose if neccessary.  Convention is to use lower diagonal
            if (x > y) {
                x = originY;
                y = originX;
                int tmp = width;
                width = height;
                height = tmp;
            }
        }

        int maxX = x + width;
        int maxY = y + height;


        if (displayOption == MainWindow.DisplayOption.PEARSON) {
            BasicMatrix bm = zd.getPearsons();
            if (bm != null) {
                ((HiCColorScale) pearsonColorScale).setMin(bm.getLowerValue());
                ((HiCColorScale) pearsonColorScale).setMax(bm.getUpperValue());
                renderMatrix(bm, originX, originY, width, height, pearsonColorScale, g);

            }
        } else {
            // Iterate through blocks overlapping visible region
            DensityFunction df = null;
            if (displayOption == MainWindow.DisplayOption.OE) {
                df = hic.getDensityFunction(zd.getZoom());
            }

            List<Block> blocks = zd.getBlocksOverlapping(x, y, maxX, maxY);

            ColorScale cs;

            if (displayOption == MainWindow.DisplayOption.OBSERVED) {
                observedColorScale = observedColorScaleMap.get(zd);
                if (observedColorScale == null) {
                    float percent90 = computePercentile(blocks, 90 );
                    observedColorScale = new ContinuousColorScale(0, percent90, Color.white, Color.red);
                    observedColorScaleMap.put(zd, observedColorScale);
                    mainWindow.updateColorSlider(0, (int) (2*percent90), (int) percent90);
                }
                cs = observedColorScale;
            } else {
                cs = oeColorScale;
            }


            for (Block b : blocks) {
                renderBlock(originX, originY, chr1, chr2, 1.0f, b, cs, df, g);
            }
        }
    }


    private float computePercentile(List<Block> blocks, double p) {

        DoubleArrayList dal = new DoubleArrayList(10000);
        for (Block b : blocks) {
            ContactRecord[] recs = b.getContactRecords();
            if (recs != null) {
                for (int i = 0; i < recs.length; i++) {
                    ContactRecord rec = recs[i];
                    if (Math.abs(rec.getBinX() - rec.getBinY()) > 1) {
                        dal.add(rec.getCounts());
                    }
                }
            }
        }


        return dal.size() == 0 ? 1 : (float) StatUtils.percentile(dal.toArray(), p);
    }


    private void renderBlock(int originX, int originY, int chr1, int chr2,
                             double colorScaleFactor,
                             Block b,
                             ColorScale colorScale, DensityFunction df, Graphics2D g) {

        MainWindow.DisplayOption displayOption = hic.getDisplayOption();
        boolean sameChr = (chr1 == chr2);

        ContactRecord[] recs = b.getContactRecords();
        if (recs != null) {
            for (int i = 0; i < recs.length; i++) {
                ContactRecord rec = recs[i];

                Color color;
                double score;
                // This is weirdly not the same as computeOE.
                if (displayOption == MainWindow.DisplayOption.OE && df != null) {
                    int x = rec.getBinX();// * binSize;
                    int y = rec.getBinY();// * binSize;
                    int dist = Math.abs(x - y);
                    double expected = df.getDensity(chr1, dist);
                    double observed = rec.getCounts();
                    // double normCounts = (rec.getCounts() / expected);
                    double normCounts = observed / expected;
                    score = normCounts;
                    score = Math.log10(score);
                } else {
                    score = rec.getCounts() / colorScaleFactor;
                }

                color = colorScale.getColor((float) score);
                int px = (rec.getBinX() - originX);
                int py = (rec.getBinY() - originY);
                g.setColor(color);
                // TODO -- need to check right bounds before drawing
                if (px > -1 && py > -1) {
                    g.fillRect(px, py, MainWindow.BIN_PIXEL_WIDTH, MainWindow.BIN_PIXEL_WIDTH);
                }

                if (sameChr && (rec.getBinX() != rec.getBinY())) {
                    px = (rec.getBinY() - originX);
                    py = (rec.getBinX() - originY);
                    if (px > -1 && py > -1) {
                        g.fillRect(px, py, MainWindow.BIN_PIXEL_WIDTH, MainWindow.BIN_PIXEL_WIDTH);
                    }
                }
            }
        }
    }

    /**
     * Used for Pearsons correlation (dense matrix).  The bitmap is drawn at 1 data point per pixel, scaling
     * happens elsewhere.
     *
     * @param originX    origin in pixels
     * @param originY    origin in pixels
     * @param rm
     * @param colorScale
     * @param g
     */
    private void renderMatrix(BasicMatrix rm, int originX, int originY, int width, int height,
                              ColorScale colorScale, Graphics g) {


        int endX = Math.min(originX + width, rm.getColumnDimension());
        int endY = Math.min(originY + height, rm.getRowDimension());

        // TODO -- need to check bounds before drawing
        for (int row = originY; row < endY; row++) {
            for (int col = originX; col < endX; col++) {

                float score = rm.getEntry(row, col);
                Color color;
                if (Float.isNaN(score)) {
                    color = Color.gray;
                } else {
                    color = score == 0 ? Color.black : colorScale.getColor((float) score);
                }
                int px = col - originX;
                int py = row - originY;
                g.setColor(color);
                g.fillRect(px, py, MainWindow.BIN_PIXEL_WIDTH, MainWindow.BIN_PIXEL_WIDTH);
                // Assuming same chromosome
                if (col != row) {
                    px = (row - originX);
                    py = (col - originY);
                    g.fillRect(px, py, MainWindow.BIN_PIXEL_WIDTH, MainWindow.BIN_PIXEL_WIDTH);
                }
            }
        }
    }


    public void setObservedRange(float min, float max) {
        if (observedColorScale == null) {
            observedColorScale = new ContinuousColorScale(min, max, Color.white, Color.red);
        }
        observedColorScale.setNegEnd(min);
        observedColorScale.setPosEnd(max);
    }
}
