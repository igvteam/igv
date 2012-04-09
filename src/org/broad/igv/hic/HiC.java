package org.broad.igv.hic;

import org.broad.igv.hic.data.*;

import javax.swing.*;

/**
 * This is the "model" class for the HiC viewer.
 *
 * @author Jim Robinson
 * @date 4/8/12
 */
public class HiC {

    MainWindow mainWindow;

    MainWindow.DisplayOption displayOption;
    Dataset dataset;
    Context xContext;
    Context yContext;
    Matrix matrix;
    MatrixZoomData zd;


    public HiC(MainWindow mainWindow) {
        this.mainWindow = mainWindow;
    }

    public MainWindow.DisplayOption getDisplayOption() {
        return displayOption;
    }

    /**
     * Change zoom level and recenter.  Triggered by the resolutionSlider, or by a double-click in the
     * heatmap panel.
     *
     * @param newZoom
     * @param centerLocationX center X location in base pairs
     * @param centerLocationY center Y location in base pairs
     */
    public void setZoom(int newZoom, final int centerLocationX, final int centerLocationY, boolean updateSlider) {

        if (newZoom < 0 || newZoom > MainWindow.MAX_ZOOM) return;

        final Chromosome chr1 = xContext.getChromosome();
        final Chromosome chr2 = yContext.getChromosome();
        final MatrixZoomData newZD = dataset.getMatrix(chr1, chr2).getObservedMatrix(newZoom);


        if (displayOption == MainWindow.DisplayOption.PEARSON && newZD.getPearsons() == null) {

            if (newZoom > 3) {
                int ans = JOptionPane.showConfirmDialog(mainWindow.resolutionSlider.getTopLevelAncestor(), "Pearson's calculation at " +
                        "this zoom will take a while.\nAre you sure you want to proceed?",
                        "Confirm calculation", JOptionPane.YES_NO_OPTION);
                if (ans == JOptionPane.NO_OPTION) {
                    mainWindow.resolutionSlider.setValue(zd.getZoom());
                    return;
                }
            }

            final DensityFunction df = mainWindow.getDensityFunction(newZoom);

            Runnable callable = new Runnable() {
                public void run() {
                    newZD.computePearsons(df);
                    mainWindow.refresh();
                    updateState2(newZD, centerLocationX, centerLocationY);
                }
            };
            mainWindow.executeLongRunningTask(callable);
            return;


        } else {
            newZD.printDescription();
            updateState2(newZD, centerLocationX, centerLocationY);
        }

    }


    private void updateState2(MatrixZoomData newZD, double centerLocationX, double centerLocationY) {

        zd = newZD;
        int newZoom = zd.getZoom();
        //showEigenvector();

        int newBinSize = zd.getBinSize();

        double xScaleMax = (double) xContext.getChrLength() / mainWindow.getHeatmapPanel().getWidth();
        double yScaleMax = (double) yContext.getChrLength() / mainWindow.getWidth();
        double scaleMax = Math.max(xScaleMax, yScaleMax);

        double scale = Math.min(newBinSize, scaleMax);

        xContext.setZoom(newZoom, scale);
        yContext.setZoom(newZoom, scale);

        center((int) centerLocationX, (int) centerLocationY);
        mainWindow.resolutionSlider.setValue(newZoom);
        mainWindow.refresh();
    }


    /**
     * Zoom to a specific rectangle.  Triggered by the alt-drag action in the heatmap panel.
     *
     * @param xBP
     * @param yBP
     * @param scale
     */
    public void zoomTo(final double xBP, final double yBP, final double scale) {

        // Find zoom level closest to prescribed scale
        int newZoom = HiCGlobals.zoomBinSizes.length - 1;
        for (int z = 1; z < HiCGlobals.zoomBinSizes.length; z++) {
            if (HiCGlobals.zoomBinSizes[z] < scale) {
                newZoom = z - 1;
                break;
            }
        }

        final Chromosome chr1 = xContext.getChromosome();
        final Chromosome chr2 = yContext.getChromosome();
        final MatrixZoomData newZD = dataset.getMatrix(chr1, chr2).getObservedMatrix(newZoom);

        if (displayOption == MainWindow.DisplayOption.PEARSON && newZD.getPearsons() == null) {
            if (newZoom > 3) {
                int ans = JOptionPane.showConfirmDialog(mainWindow.resolutionSlider.getTopLevelAncestor(), "Pearson's calculation at " +
                        "this zoom will take a while.\nAre you sure you want to proceed?",
                        "Confirm calculation", JOptionPane.YES_NO_OPTION);
                if (ans == JOptionPane.NO_OPTION) {
                    return;
                }
            }

            final DensityFunction df = mainWindow.getDensityFunction(newZoom);
            Runnable callable = new Runnable() {
                public void run() {
                    newZD.computePearsons(df);
                    mainWindow.refresh();
                    updateState(newZD, scale, xBP, yBP);
                }
            };
            mainWindow.executeLongRunningTask(callable);
            return;


        } else {
            updateState(newZD, scale, xBP, yBP);

        }
    }


    private void updateState(MatrixZoomData newZD, double scale, double xBP, double yBP) {
        zd = newZD;
        xContext.setZoom(zd.getZoom(), scale);
        yContext.setZoom(zd.getZoom(), scale);
        xContext.setOrigin((int) xBP);
        yContext.setOrigin((int) yBP);
        mainWindow.resolutionSlider.setValue(zd.getZoom());
        mainWindow.refresh();
    }

    public void center(int centerLocationX, int centerLocationY) {
        double w = (mainWindow.getHeatmapPanel().getWidth() * xContext.getScale());
        int newX = (int) (centerLocationX - w / 2);
        double h = (mainWindow.getHeatmapPanel().getHeight() * yContext.getScale());
        int newY = (int) (centerLocationY - h / 2);
        moveTo(newX, newY);
    }


    public void moveBy(int dx, int dy) {
        final int newX = xContext.getOrigin() + dx;
        final int newY = yContext.getOrigin() + dy;
        moveTo(newX, newY);
    }

    private void moveTo(int newX, int newY) {
        final double bpWidthX = xContext.getScale() * mainWindow.getHeatmapPanel().getWidth();
        int maxX = (int) (xContext.getChrLength() - bpWidthX);
        final double bpWidthY = yContext.getScale() * mainWindow.getHeatmapPanel().getHeight();
        int maxY = (int) (yContext.getChrLength() - bpWidthY);

        int x = Math.max(0, Math.min(maxX, newX));
        int y = Math.max(0, Math.min(maxY, newY));

        xContext.setOrigin(x);
        yContext.setOrigin(y);

//        String locus1 = "chr" + (xContext.getChromosome().getName()) + ":" + x + "-" + (int) (x + bpWidthX);
//        String locus2 = "chr" + (yContext.getChromosome().getName()) + ":" + x + "-" + (int) (y + bpWidthY);
//        IGVUtils.sendToIGV(locus1, locus2);

        mainWindow.repaint();
    }


    public void setDisplayOption(String strOption) {

        MainWindow.DisplayOption newDisplay = MainWindow.DisplayOption.getByValue(strOption);

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
                final DensityFunction df = mainWindow.getDensityFunction(zd.getZoom());
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

    public void reset() {
        matrix = null;
        zd = null;
    }
}
