package org.igv.sashimi;

import org.igv.event.IGVEventBus;
import org.igv.feature.Locus;
import org.igv.ui.panel.ReferenceFrame;

/**
 * Horizontal coordinates of a Sashimi plot.  Pan and zoom act on the plot frame, whose origin and scale are in plot
 * coordinates, i.e. genomic coordinates transformed by a {@link SashimiCoordinateMap}.  The data frame spans the same
 * view in genomic coordinates, and is used to load data and to find features under the mouse.
 */
public class SashimiView {

    private final ReferenceFrame plotFrame;
    private final ReferenceFrame dataFrame;
    private SashimiCoordinateMap map = SashimiCoordinateMap.identity();

    public SashimiView(ReferenceFrame frame, IGVEventBus eventBus) {
        this.plotFrame = new ReferenceFrame(frame, eventBus);
        this.dataFrame = new DataFrame(frame);
    }

    public ReferenceFrame getPlotFrame() {
        return plotFrame;
    }

    public ReferenceFrame getDataFrame() {
        return dataFrame;
    }

    public SashimiCoordinateMap getMap() {
        return map;
    }

    /**
     * Replace the coordinate map, keeping the genomic range in view.  The plot frame posts a view change event.
     */
    public void setMap(SashimiCoordinateMap newMap) {
        double start = toGenomic(0);
        double end = toGenomic(plotFrame.getWidthInPixels());
        map = newMap;
        showGenomicRange(start, end);
    }

    /**
     * Re-fit the view to a new pixel width, keeping the genomic range in view.  Called when the window is resized;
     * without it the plot is drawn at the scale computed for the old width.
     */
    public void setWidthInPixels(int widthInPixels) {
        if (widthInPixels <= 0 || widthInPixels == plotFrame.getWidthInPixels()) {
            return;
        }
        double start = toGenomic(0);
        double end = toGenomic(plotFrame.getWidthInPixels());
        plotFrame.setWidthInPixels(widthInPixels);
        dataFrame.setWidthInPixels(widthInPixels);
        showGenomicRange(start, end);
    }

    /**
     * Fit the given genomic range to the plot frame's width.  The plot frame posts a view change event.
     */
    private void showGenomicRange(double start, double end) {
        int plotStart = (int) Math.round(map.toPlot(start));
        // Keep a non-zero width -- a heavily shrunk view can round to a single coordinate
        int plotEnd = Math.max(plotStart + 1, (int) Math.round(map.toPlot(end)));
        plotFrame.jumpTo(new Locus(plotFrame.getChrName(), plotStart, plotEnd));
        syncDataFrame();
    }

    /**
     * Set the data frame to the genomic range of the plot frame's view.
     */
    public void syncDataFrame() {
        int start = (int) Math.floor(toGenomic(0));
        int end = (int) Math.ceil(toGenomic(plotFrame.getWidthInPixels()));
        dataFrame.jumpTo(new Locus(plotFrame.getChrName(), start, end));
    }

    public double toPixel(double position) {
        return (map.toPlot(position) - plotFrame.getOrigin()) / plotFrame.getScale();
    }

    public double toGenomic(double pixel) {
        return map.toGenomic(plotFrame.getOrigin() + pixel * plotFrame.getScale());
    }

    /**
     * Genomic frame whose screen positions follow the plot frame through the coordinate map.
     */
    private class DataFrame extends ReferenceFrame {

        DataFrame(ReferenceFrame frame) {
            super(frame, new IGVEventBus());
        }

        @Override
        public double getChromosomePosition(int screenPosition) {
            return toGenomic(screenPosition);
        }
    }
}
