/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.Interval;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ObservableForObject;

import java.util.Observer;


/**
 * @author jrobinso
 */
public class ReferenceFrame {

    private static Logger log = Logger.getLogger(ReferenceFrame.class);

    /**
     * The nominal viewport width in pixels.
     */
    public static int binsPerTile = 700;

    private String name;

    int pixelX;

    protected int widthInPixels;


    /**
     * The chromosome currently in view
     */
    protected String chrName = "chrAll";

    /**
     * The minimum zoom level for the current screen size + chromosome combination.
     */
    private int minZoom = 0;

    /**
     * The current zoom level.  Zoom level -1 corresponds to the whole
     * genome view (chromosome "all")
     */
    protected int zoom = minZoom;
    /**
     * The maximum zoom level.  Set to prevent integer overflow.  This is a function
     * of chromosom length.
     */
    public static int maxZoom = 23;
    /**
     * The number of tiles for this zoom level, = 2^zoom
     */
    private int nTiles = 1;
    /**
     * The maximum virtual pixel value.
     */
    private double maxPixel;
    /**
     * The origin in bp
     */
    protected double origin = 0;

    /**
     * The location (x axis) locationScale in base pairs / virtual pixel
     */
    protected double locationScale;

    /**
     * Minimum allowed scale in base pairs / pixel
     */
    private double minScale = 1.0 / 14;

    protected boolean locationScaleValid = false;
    protected Locus initialLocus;
    protected int setEnd = 0;

    protected ObservableForObject<String> chromoObservable;


    public ReferenceFrame(String name) {
        this.name = name;
        Genome genome = getGenome();
        this.chrName = genome == null ? "" : genome.getHomeChromosome();
        chromoObservable = new ObservableForObject<String>(chrName);
    }


    public ReferenceFrame(ReferenceFrame otherFrame) {
        this.chrName = otherFrame.chrName;
        this.initialLocus = otherFrame.initialLocus;
        this.locationScale = otherFrame.locationScale;
        this.locationScaleValid = otherFrame.locationScaleValid;
        this.maxPixel = otherFrame.maxPixel;
        this.minScale = otherFrame.minScale;
        this.minZoom = otherFrame.minZoom;
        this.name = otherFrame.name;
        this.nTiles = otherFrame.nTiles;
        this.origin = otherFrame.origin;
        this.pixelX = otherFrame.pixelX;
        this.setEnd = otherFrame.setEnd;
        this.widthInPixels = otherFrame.widthInPixels;
        this.zoom = otherFrame.zoom;
        chromoObservable = new ObservableForObject<String>(chrName);
    }


    public synchronized void setBounds(int x, int w) {
        this.pixelX = x;
        if (w != widthInPixels) {
            widthInPixels = w;
            computeLocationScale();
        }
    }

    public int getMidpoint() {
        return pixelX + widthInPixels / 2;
    }

    /**
     * Compute the maximum zoom level, which is a function of chromosome length.
     */
    public void computeMaxZoom() {
        Genome genome = getGenome();
        // Compute max zoom.  Assume window size @ max zoom of ~ 50 bp
        if (genome != null && chrName != null && genome.getChromosome(chrName) != null) {
            if (chrName.equals(Globals.CHR_ALL)) {
                maxZoom = 0;
            } else {
                int chrLength = genome.getChromosome(chrName).getLength();

                maxZoom = (int) (Math.log(chrLength / 50.0) / Globals.log2) + 1;

            }

            if (zoom > maxZoom) {
                zoomTo(maxZoom);
            }
        }

    }

    /**
     * Change the zoom level, keeping the center the same
     * zoom events "release" the frame, enabling pan and zoom
     *
     * @param newZoom
     */
    private void zoomTo(int newZoom) {

        zoom = Math.min(maxZoom, newZoom);
        nTiles = (int) Math.pow(2, Math.max(minZoom, zoom));
        maxPixel = getTilesTimesBinsPerTile();
        invalidateLocationScale();

        // TODO -- do this with events,
        IGV.repaintPanelsHeadlessSafe();
    }

    /**
     * This setter is provided to restore state, it does not do a "zoom" action
     *
     * @param z
     */
    public void setZoom(int z) {
        if (z > 0)
            this.zoom = z;
    }

    public void incrementZoom(int increment) {
        double currentCenter = getGenomeCenterPosition();
        zoomBy(increment, currentCenter);
    }

    public void zoomAndCenterAdjusted(int newZoom) {
        double currentCenter = getGenomeCenterPosition();
        zoomTo(newZoom, currentCenter);
    }

    protected double getGenomeCenterPosition() {
        return origin + ((widthInPixels / 2) * getScale());
    }

    public synchronized void zoomBy(final int zoomFactor, final double newCenter) {

        if (FrameManager.isGeneListMode()) {
            double f = Math.pow(2.0, zoomFactor);
            setLocationScale(Math.max(minScale, getScale() / f));
            double newOrigin = Math.round(newCenter - ((widthInPixels / 2) * getScale()));
            setOrigin(newOrigin);
            double end = setEnd > 0 ? setEnd : getEnd();
            imputeZoom(origin, end);
            setEnd = -1;
        } else {
            int newZoom = Math.max(0, zoom + zoomFactor);
            zoomTo(newZoom, newCenter);
        }
    }

    public synchronized void zoomTo(final int newZoom, final double newCenter) {

        if (chrName.equals(Globals.CHR_ALL)) {
            chrName = getGenome().getHomeChromosome();
        }

        if (chrName.equals(Globals.CHR_ALL)) {
            // Translate the location to chromosome number
            jumpToChromosomeForGenomeLocation(newCenter);
            IGV.getInstance().chromosomeChangeEvent(chrName);
        } else {
            if (zoom != newZoom) {
                zoomTo(newZoom);
                double newLocationScale = getScale();
                // Adjust origin so newCenter is centered
                double newOrigin = Math.round(newCenter - ((widthInPixels / 2) * newLocationScale));
                setOrigin(newOrigin);

            }
        }

        recordHistory();

        // This is a hack,  this is not a drag event but is a "jump"
        DragEventManager.getInstance().dragStopped();

    }

    public void recordHistory() {
        IGV.getInstance().getSession().getHistory().push(getFormattedLocusString(), zoom);
    }

    private void jumpToChromosomeForGenomeLocation(double locationMB) {
        double startMB = 0;

        for (String chr : getGenome().getLongChromosomeNames()) {
            double endMB = startMB + getGenome().getChromosome(chr).getLength() / 1000.0;

            if ((locationMB > startMB) && (locationMB <= endMB)) {
                this.setChromosomeName(chr);
                break;
            }

            startMB = endMB;
        }
    }

    public void shiftOriginPixels(int delta) {
        double shiftBP = delta * getScale();
        setOrigin(origin + shiftBP);
    }

    public void snapToGrid() {
        setOrigin(Math.round(origin));
    }

    public void centerOnLocation(String chr, double chrLocation) {
        if (!chrName.equals(chr)) {
            chrName = chr;
            computeMaxZoom();
            if (zoom > maxZoom) {
                zoomTo(maxZoom);
            }
        }
        centerOnLocation(chrLocation);
    }

    public void centerOnLocation(double chrLocation) {
        double windowWidth = (widthInPixels * getScale()) / 2;
        setOrigin(Math.round(chrLocation - windowWidth));
        recordHistory();
    }

    public boolean windowAtEnd() {
        double windowLengthBP = widthInPixels * getScale();
        return origin + windowLengthBP + 1 > getChromosomeLength();

    }

    public void setOrigin(double position) {
        setOrigin(position, true);

    }

    /**
     * Set the origin of the frame, guarding against chromosome boundaries
     *
     * @param position
     * @param repaint
     */
    public void setOrigin(double position, boolean repaint) {

        int windowLengthBP = (int) (widthInPixels * getScale());
        double newOrigin;
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_SOFT_CLIPPED)) {
            newOrigin = Math.max(-1000, Math.min(position, getChromosomeLength() + 1000 - windowLengthBP));
        } else {
            newOrigin = Math.max(0, Math.min(position, getChromosomeLength() - windowLengthBP));
        }
        //double delta = newOrigin - origin;
        origin = newOrigin;

        if (repaint) {
            IGV.repaintPanelsHeadlessSafe();
        }
    }

    /**
     * Move the frame to the specified position. New zoom is calculated
     * based on limits.
     *
     * @param chr
     * @param start
     * @param end
     */
    public void jumpTo(String chr, int start, int end) {
        Locus locus = new Locus(chr, start, end);
        this.jumpTo(locus);
    }

    public void jumpTo(Locus locus) {
        String chr = locus.getChr();
        int start = locus.getStart();
        int end = locus.getEnd();

        Genome genome = getGenome();
        if (chr != null) {
            if (genome.getChromosome(chr) == null && !chr.contains(Globals.CHR_ALL)) {
                MessageUtils.showMessage(chr + " is not a valid chromosome.");
                return;
            }
        }

        this.initialLocus = locus;
        Chromosome chromosome = genome == null ? null : genome.getChromosome(chr);
        if (chromosome != null) {
            end = Math.min(chromosome.getLength(), end);
        }

        synchronized (this) {
            this.chrName = chr;
            if (start >= 0) {
                imputeZoom(start, end);
                if (widthInPixels > 0) {
                    setLocationScale(((double) (end - start)) / widthInPixels);
                }
                this.setEnd = locus.getEnd();
                origin = start;
            }
        }

        if (log.isDebugEnabled()) {
            log.debug("Data panel width = " + widthInPixels);
            log.debug("New start = " + (int) origin);
            log.debug("New end = " + (int) getEnd());
            log.debug("New center = " + (int) getCenter());
            log.debug("Scale = " + locationScale);
        }

        //Mostly for testing
        IGV.repaintPanelsHeadlessSafe();
    }

    protected void imputeZoom(double start, double end) {
        int z = (int) (Math.log(getChromosomeLength() / (end - start)) / Globals.log2) + 1;
        if (z != this.zoom) {
            zoom = Math.min(maxZoom, Math.max(minZoom, z));
            nTiles = (int) Math.pow(2, zoom);
            maxPixel = getTilesTimesBinsPerTile();
        }
        if (IGV.hasInstance())
            IGV.getInstance().repaintStatusAndZoomSlider();
    }

    protected Genome getGenome() {
        return GenomeManager.getInstance().getCurrentGenome();
    }

    public double getOrigin() {
        return origin;
    }

    public double getCenter() {
        return origin + getScale() * widthInPixels / 2;
    }

    public double getEnd() {
        return origin + getScale() * widthInPixels;
    }

    public int getZoom() {
        return zoom;
    }

    /**
     * Return the maximum zoom level
     *
     * @return
     */
    public int getMaxZoom() {
        return maxZoom;
    }

    public int getAdjustedZoom() {
        return zoom - minZoom;
    }

    public double getMaxPixel() {
        return maxPixel;
    }

    /**
     * Change the frame to the specified chromosome.
     *
     * @param name
     * @param force
     */
    public void setChromosomeName(String name, boolean force) {

        if ((chrName == null) || !name.equals(chrName) || force) {
            chrName = name;
            origin = 0;
            setEnd = -1;
            zoomTo(0);
            computeMaxZoom();

            chromoObservable.setChangedAndNotify();
        }
    }

    public void addObserver(Observer observer) {
        chromoObservable.addObserver(observer);
    }

    public void deleteObserver(Observer observer) {
        chromoObservable.deleteObserver(observer);
    }

    public void deleteObservers() {
        chromoObservable.deleteObservers();
    }

    /**
     * @param name
     */
    public void setChromosomeName(String name) {
        setChromosomeName(name, false);
    }


    public String getChrName() {
        return chrName;
    }


    // TODO -- this parameter shouldn't be stored here.  Maybe in a specialized
    // layout manager?

    /**
     * Width in pixels
     */

    public int getWidthInPixels() {
        return widthInPixels;
    }

    /**
     * Return the current locationScale in base pairs / pixel
     *
     * @return
     */
    public synchronized double getScale() {
        if ((locationScale <= 0) || !locationScaleValid) {
            computeLocationScale();
        }

        return locationScale;
    }


    public void invalidateLocationScale() {
        //Thread.dumpStack();
        locationScaleValid = false;
    }

    private synchronized void computeLocationScale() {
        Genome genome = getGenome();

        if (genome != null) {

            if (setEnd > 0 && widthInPixels > 0) {
                setLocationScale((setEnd - origin) / widthInPixels);
                imputeZoom(origin, setEnd);
                setEnd = -1;
            } else {
                computeMinZoom();
                double virtualPixelSize = getTilesTimesBinsPerTile();
                double nPixel = Math.max(virtualPixelSize, widthInPixels);
                setLocationScale(((double) getChromosomeLength()) / nPixel);
            }
        }
    }

    /**
     * Compute the minimum zoom level for the data panel width.  This is defined as the maximum
     * zoom level for which all data bins will fit in the window without loss of
     * data,  i.e. the maximum zoom level for which nBins < nPixels.  The number
     * of bins is defined as
     * nBins =  2^z
     * so minZoom is the value z such that nBins < dataPanelWidth
     */
    protected void computeMinZoom() {
        if (this.chrName.equals(Globals.CHR_ALL)) {
            minZoom = 0;
        } else {
            minZoom = Math.max(0, (int) (Math.log((widthInPixels / binsPerTile)) / Globals.log2));

            if (zoom < minZoom) {
                zoom = minZoom;
                nTiles = (int) Math.pow(2, zoom);
                maxPixel = getTilesTimesBinsPerTile();
            }
        }

    }

    /**
     * Return the chromosome position corresponding to the pixel index.  The
     * pixel index is the pixel "position" translated by -origin.
     *
     * @param screenPosition
     * @return
     */
    public double getChromosomePosition(int screenPosition) {
        return origin + getScale() * screenPosition;
    }

    /**
     * Return the screen position corresponding to the chromosomal position.
     *
     * @param chromosomePosition
     * @return
     */
    public int getScreenPosition(double chromosomePosition) {
        return (int) ((chromosomePosition - origin) / getScale());
    }


    public Chromosome getChromosome() {
        Genome genome = getGenome();
        if (genome == null) {
            return null;
        }

        return genome.getChromosome(chrName);
    }


    public int getChromosomeLength() {
        Genome genome = getGenome();

        if (genome == null) {
            return 1;
        }

        if (chrName.equals("All")) {

            // TODO -- remove the hardcoded unit divider ("1000")
            return (int) (genome.getNominalLength() / 1000);

            // return genome.getLength();
        } else {
            if (getChromosome() == null) {
                log.error("Null chromosome: " + chrName);
                if (genome == null || genome.getChromosomes().size() == 0) {
                    return 1;
                } else {
                    return genome.getChromosomes().iterator().next().getLength();
                }
            }

            return getChromosome().getLength();
        }
    }


    public double getTilesTimesBinsPerTile() {
        return (double) nTiles * (double) binsPerTile;
    }

    /**
     * Get the UCSC style locus string corresponding to the current view.  THe UCSC
     * conventions are followed for coordinates,  specifically the internal representation
     * is "zero" based (first base is numbered 0) but the display representation is
     * "one" based (first base is numbered 1).   Consequently 1 is added to the
     * computed positions.
     *
     * @return
     */
    public String getFormattedLocusString() {

        if (zoom == 0) {
            return getChrName();
        } else {

            Range range = getCurrentRange();
            return Locus.getFormattedLocusString(range.getChr(), range.getStart(), range.getEnd());
        }
    }

    public Range getCurrentRange() {
        int start = 0;
        int end = widthInPixels;
        int startLoc = (int) getChromosomePosition(start) + 1;
        int endLoc = (int) getChromosomePosition(end);
        Range range = new Range(getChrName(), startLoc, endLoc);
        return range;
    }

    protected synchronized void setLocationScale(double locationScale) {
        if (log.isDebugEnabled()) {
            log.debug("Set location scale: " + locationScale + "  origin=" + origin);
        }
        this.locationScale = locationScale;
        locationScaleValid = true;

    }

    public void reset() {
        jumpTo(FrameManager.getLocus(name));
    }

    public String getName() {
        return name;
    }

    public Locus getInitialLocus() {
        return initialLocus;
    }


    public int getMinZoom() {
        return minZoom;
    }

    public boolean isExomeMode() {
        return false;
    }

    public void setName(String name) {
        this.name = name;
    }

    public boolean overlaps(Interval interval) {
        return this.getChrName().equals(interval.getChr()) && this.getOrigin() <= interval.getEnd() && this.getEnd() >= interval.getStart();
    }

    public static class Range {
        private String chr;
        private int start;
        private int end;

        public Range(String chr, int start, int end) {
            this.chr = chr;
            this.start = start;
            this.end = end;
        }

        public String getChr() {
            return chr;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public int getLength() {
            return end - start;
        }
    }

}
