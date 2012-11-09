package org.broad.igv.hic;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.data.*;

import javax.swing.*;
import java.awt.*;
import java.util.Map;

/**
 * This is the "model" class for the HiC viewer.
 *
 * @author Jim Robinson
 * @date 4/8/12
 */
public class HiC {


    private int cursorX;

    public Unit getUnit() {
        return unit;
    }

    public void setCursorX(int x) {
        this.cursorX = x;

    }

    public int getCursorX() {
        return cursorX;
    }

    public enum Unit {BP, FRAG}

    ;

    MainWindow mainWindow;
    Unit unit = Unit.BP;
    MainWindow.DisplayOption displayOption;
    Dataset dataset;
    private Chromosome[] chromosomes;

    public Context xContext;
    Context yContext;

    Matrix matrix;
    public MatrixZoomData zd;


    public HiC(MainWindow mainWindow) {
        this.mainWindow = mainWindow;
    }

    public MainWindow.DisplayOption getDisplayOption() {
        return displayOption;
    }


    public boolean isWholeGenome() {
        return xContext != null && xContext.getChromosome().getName().equals("All");
    }


    public void reset() {
        matrix = null;
        zd = null;
    }


    public void setChromosomes(Chromosome[] chromosomes) {
        this.chromosomes = chromosomes;
    }

    public Chromosome[] getChromosomes() {
        return chromosomes;
    }


    /**
     * Return the expected density function for the current zoom level
     *
     * @return
     */
    public DensityFunction getDensityFunction() {

        if (zd == null) return null;

        return getDensityFunction(zd.getZoom());
    }

    /**
     * Return the expected density function for a specified zoom level.
     *
     * @param zoom
     * @return density function, or null if expected densities are not available
     */
    public DensityFunction getDensityFunction(int zoom) {

        return dataset.getDensityFunction(zoom);
    }


    public void setUnit(Unit unit) {
        this.unit = unit;
    }

    /**
     * Change zoom level and recenter.  Triggered by the resolutionSlider, or by a double-click in the
     * heatmap panel.
     *
     * @param newZoom
     * @param genomePositionX center X in base pairs
     * @param genomePositionY center Y in base pairs
     */
    public void setZoom(int newZoom, final double genomePositionX, final double genomePositionY) {

        if (newZoom < 0 || newZoom > dataset.getNumberZooms(unit)) return;

        final Chromosome chr1 = xContext.getChromosome();
        final Chromosome chr2 = yContext.getChromosome();


        int zoomDataIdx = unit == HiC.Unit.BP ? newZoom : newZoom + 9;

        final MatrixZoomData newZD = dataset.getMatrix(chr1, chr2).getObservedMatrix(zoomDataIdx);
        if (newZD == null) {
            JOptionPane.showMessageDialog(mainWindow, "Sorry, this zoom is not available", "Zoom unavailable", JOptionPane.WARNING_MESSAGE);
            return;
        }

        if (displayOption == MainWindow.DisplayOption.PEARSON && newZD.getPearsons() == null) {

            if (newZoom > 3) {
                int ans = JOptionPane.showConfirmDialog(mainWindow, "Pearson's calculation at " +
                        "this resolution will take a while.\nAre you sure you want to proceed?",
                        "Confirm calculation", JOptionPane.YES_NO_OPTION);
                if (ans == JOptionPane.NO_OPTION) {
                    mainWindow.updateZoom(unit, zd.getZoom());
                    return;
                }
            }

            final DensityFunction df = getDensityFunction(newZoom);
            Runnable callable = new Runnable() {
                public void run() {
                    newZD.computePearsons(df);
                    //mainWindow.updateEigenvectorTrack();
                    mainWindow.refresh();
                    updateState2(newZD, genomePositionX, genomePositionY);
                }
            };
            mainWindow.executeLongRunningTask(callable);
        } else {
            //newZD.printDescription();
            updateState2(newZD, genomePositionX, genomePositionY);
        }

    }


    private void updateState2(MatrixZoomData newZD, double genomePositionX, final double genomePositionY) {

        zd = newZD;
        int newZoom = zd.getZoom();
        //showEigenvector();

        //int newBinSize = zd.getBinSize();

        // double xScaleMax = (double) xContext.getChrLength() / mainWindow.getHeatmapPanel().getWidth();
        // double yScaleMax = (double) yContext.getChrLength() / mainWindow.getWidth();
        // double scaleMax = Math.max(xScaleMax, yScaleMax);

        //  double scale = Math.min(newBinSize, scaleMax);

        xContext.setZoom(newZoom);
        yContext.setZoom(newZoom);

        int xBinCount = zd.getxGridAxis().getBinCount();
        int yBinCount = zd.getyGridAxis().getBinCount();
        int maxBinCount = Math.max(xBinCount, yBinCount);
        double scalefactor = Math.max(1.0, (double) mainWindow.getHeatmapPanel().getWidth() / maxBinCount);
        xContext.setScaleFactor(scalefactor);
        yContext.setScaleFactor(scalefactor);

        //Point binPosition = zd.getBinPosition(genomePositionX, genomePositionY);
        int binX = zd.getxGridAxis().getBinNumberForGenomicPosition((int) genomePositionX);
        int binY = zd.getyGridAxis().getBinNumberForGenomicPosition((int) genomePositionY);

        center(binX, binY);
        mainWindow.updateZoom(unit, newZoom);

        mainWindow.refresh();
    }


    /**
     * Zoom to a specific rectangle.  Triggered by the alt-drag action in the heatmap panel.
     */
    public void zoomTo(final int xBP0, final int yBP0, final int xBP1, int yBP1, double genomicScale) {

        // Find zoom level closest to prescribed scale
        int newZoom = dataset.getNumberZooms(unit) - 1;
        for (int z = 1; z < dataset.getNumberZooms(unit); z++) {
            if (dataset.getZoom(unit, z) < genomicScale) {
                newZoom = z - 1;
                break;
            }
        }

        final Chromosome chr1 = xContext.getChromosome();
        final Chromosome chr2 = yContext.getChromosome();


        int zoomDataIdx = unit == HiC.Unit.BP ? newZoom : newZoom + 9;

        final MatrixZoomData newZD = dataset.getMatrix(chr1, chr2).getObservedMatrix(zoomDataIdx);


        int binX0 = newZD.getxGridAxis().getBinNumberForGenomicPosition((int) xBP0);
        int binY0 = newZD.getyGridAxis().getBinNumberForGenomicPosition((int) yBP0);
        final int binXMax = newZD.getxGridAxis().getBinNumberForGenomicPosition((int) xBP1);
        final int binYMax = newZD.getyGridAxis().getBinNumberForGenomicPosition((int) yBP1);

        final int binX = binX0;
        final int binY = binY0;
        final double xScale = ((double) mainWindow.getHeatmapPanel().getWidth()) / (binXMax - binX);
        final double yScale = ((double) mainWindow.getHeatmapPanel().getHeight()) / (binYMax - binY);
        final double scaleFactor = Math.max(1, Math.min(xScale, yScale));

        if (displayOption == MainWindow.DisplayOption.PEARSON && newZD.getPearsons() == null) {
            if (newZoom > 3) {
                int ans = JOptionPane.showConfirmDialog(mainWindow, "Pearson's calculation at " +
                        "this zoom will take a while.\nAre you sure you want to proceed?",
                        "Confirm calculation", JOptionPane.YES_NO_OPTION);
                if (ans == JOptionPane.NO_OPTION) {
                    return;
                }
            }

            final DensityFunction df = getDensityFunction(newZoom);
            Runnable callable = new Runnable() {
                public void run() {
                    newZD.computePearsons(df);
                    updateState(newZD, binX, binY, scaleFactor);
                }
            };
            mainWindow.executeLongRunningTask(callable);
        } else {
            updateState(newZD, binX, binY, scaleFactor);

        }
    }


    private void updateState(MatrixZoomData newZD, int binX, int binY, double scaleFactor) {
        zd = newZD;

        xContext.setZoom(zd.getZoom(), scaleFactor);
        yContext.setZoom(zd.getZoom(), scaleFactor);

        xContext.setBinOrigin(binX);
        yContext.setBinOrigin(binY);

        mainWindow.updateZoom(unit, zd.getZoom());
        mainWindow.refresh();
    }

    public void center(int binX, int binY) {
        double w = mainWindow.getHeatmapPanel().getWidth();
        int newX = (int) (binX - w / 2);
        double h = mainWindow.getHeatmapPanel().getHeight();
        int newY = (int) (binY - h / 2);
        moveTo(newX, newY);
    }


    public void moveBy(int dx, int dy) {
        final int newX = xContext.getBinOrigin() + dx;
        final int newY = yContext.getBinOrigin() + dy;
        moveTo(newX, newY);
    }

    private void moveTo(int newBinX, int newBinY) {

        final int w = mainWindow.getHeatmapPanel().getWidth();
        int maxX = zd.getxGridAxis().getBinCount() - w;

        final int h = mainWindow.getHeatmapPanel().getHeight();
        int maxY = zd.getyGridAxis().getBinCount() - h;

        int x = Math.max(0, Math.min(maxX, newBinX));
        int y = Math.max(0, Math.min(maxY, newBinY));

        xContext.setBinOrigin(x);
        yContext.setBinOrigin(y);

//        String locus1 = "chr" + (xContext.getChromosome().getName()) + ":" + x + "-" + (int) (x + bpWidthX);
//        String locus2 = "chr" + (yContext.getChromosome().getName()) + ":" + x + "-" + (int) (y + bpWidthY);
//        IGVUtils.sendToIGV(locus1, locus2);

        mainWindow.repaint();
    }


    public void setDisplayOption(MainWindow.DisplayOption newDisplay) {

        if (this.displayOption != newDisplay) {

            if (newDisplay == MainWindow.DisplayOption.PEARSON && zd.getPearsons() == null) {
                if (zd.getZoom() > 3) {
                    int ans = JOptionPane.showConfirmDialog(mainWindow, "Pearson's calculation at " +
                            "this zoom will take a while.\nAre you sure you want to proceed?",
                            "Confirm calculation", JOptionPane.YES_NO_OPTION);
                    if (ans == JOptionPane.NO_OPTION) {
                        return;
                    }
                }

                this.displayOption = newDisplay;
                final DensityFunction df = getDensityFunction(zd.getZoom());
                Runnable callable = new Runnable() {
                    public void run() {
                        zd.computePearsons(df);
                        mainWindow.refresh();
                    }
                };
                mainWindow.executeLongRunningTask(callable);
                return;


            }
            this.displayOption = newDisplay;
            mainWindow.refresh();

        }
    }

    public double[] getEigenvector(final int n) {
        if (zd == null) return null;
        double[] eigenvector = zd.getEigenvector();
        if (eigenvector == null) {
            final DensityFunction df = getDensityFunction();
            if (df != null) {
                if (zd.getZoom() > 2) {
                    String str = "Eigenvector calculation requires Pearson's correlation matrix.\n";
                    str += "At this zoom, calculation might take a while.\n";
                    str += "Are you sure you want to proceed?";
                    int ans = JOptionPane.showConfirmDialog(mainWindow, str, "Confirm calculation", JOptionPane.YES_NO_OPTION);
                    if (ans == JOptionPane.NO_OPTION) {
                        // mainWindow.setViewEigenvector(false);
                        return null;
                    }
                }
                Runnable runnable = new Runnable() {
                    public void run() {
                        zd.computeEigenvector(df, n);
                        mainWindow.updateEigenvectorTrack();
                    }
                };
                mainWindow.executeLongRunningTask(runnable);

            }
        }

        return eigenvector;
    }


}
