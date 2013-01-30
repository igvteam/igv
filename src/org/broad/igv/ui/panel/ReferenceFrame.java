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

import com.google.common.eventbus.AsyncEventBus;
import com.google.common.eventbus.EventBus;
import com.google.common.eventbus.Subscribe;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.Interval;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.event.ZoomChange;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;
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

    /**
     * The chromosome currently in view
     */
    protected String chrName = "chrAll";

    /**
     * The minimum zoom level for the current screen size + chromosome combination.
     */
    private int minZoom = 0;

    /**
     * The maximum zoom level.  Set to prevent integer overflow.  This is a function
     * of chromosome length.
     */
    public int maxZoom = 23;

    /**
     * Minimum allowed scale in base pairs / pixel
     * TODO Does this ever change, couldn't it be static and/or final?
     */
    private double minScale = 1.0 / 14;

    /**
     * The current zoom level.  Zoom level -1 corresponds to the whole
     * genome view (chromosome "all")
     */
    protected int zoom = minZoom;


    /**
     * X location of the frame in pixels
     */
    volatile int pixelX;

    /**
     * Width of the frame in pixels
     */
    protected int widthInPixels;

    /**
     * The number of tiles for this zoom level, = 2^zoom
     */
    private double nTiles = 1;

    /**
     * The maximum virtual pixel value.
     */
    //private double maxPixel;

    /**
     * The origin in bp
     */
    protected volatile double origin = 0;

    /**
     * The location (x axis) locationScale in base pairs / virtual pixel
     */
    protected volatile double locationScale;

    protected Locus initialLocus;

    /**
     * A temporary holder. We set the end location and then
     * later zoom/origin/etc. calculations are made
     */
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

    private EventBus eventBus;
    public EventBus getEventBus(){
        if(eventBus == null){
            eventBus = new AsyncEventBus(LongRunningTask.getThreadExecutor());
            eventBus.register(this);
        }
        //if(eventBus == null) eventBus = new EventBus("IGV");
        return eventBus;
    }
    /**
     * Set the position and width of the frame, in pixels
     * @param pixelX
     * @param widthInPixels
     */
    public synchronized void setBounds(int pixelX, int widthInPixels) {
        this.pixelX = pixelX;
        if (this.widthInPixels != widthInPixels) {
            this.widthInPixels = widthInPixels;
            computeLocationScale();
            computeZoom();
        }
    }

    /**
     * Sets zoom level and recomputes scale, iff newZoom != oldZoom
     * min/maxZoom are recalculated and respected,
     * and the locationScale is recomputed
     * @param newZoom
     */
    protected void setZoom(int newZoom){
        if(zoom != newZoom){
            synchronized (this){
                setZoomWithinLimits(newZoom);
                computeLocationScale();
            }
        }
    }


    /**
     * Set the origin of the frame, guarding against chromosome boundaries
     *
     * @param position
     */
    public void setOrigin(double position) {
        int windowLengthBP = (int) (widthInPixels * getScale());
        double newOrigin;
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_SOFT_CLIPPED)) {
            newOrigin = Math.max(-1000, Math.min(position, getChromosomeLength() + 1000 - windowLengthBP));
        } else {
            newOrigin = Math.max(0, Math.min(position, getChromosomeLength() - windowLengthBP));
        }
        origin = newOrigin;
    }


    private synchronized void setZoomWithinLimits(int newZoom){
        computeMinZoom();
        computeMaxZoom();
        zoom = Math.max(minZoom, Math.min(maxZoom, newZoom));
        nTiles = Math.pow(2, zoom);
    }

    /**
     * Increment the zoom level by {@code zoomIncrement}, leaving
     * the center the same
     * @param zoomIncrement
     */
    public void doZoomIncrement(int zoomIncrement) {
        double currentCenter = getGenomeCenterPosition();
        doIncrementZoom(zoomIncrement, currentCenter);
    }

    /**
     * Set the zoom level to {@code newZoom}, leaving
     * the center the same
     * @param newZoom
     */
    public void doSetZoom(int newZoom) {
        double currentCenter = getGenomeCenterPosition();
        doSetZoomCenter(newZoom, currentCenter);
    }

    @Subscribe
    public void receiveZoomChange(ZoomChange.Cause e) {
        doSetZoom(e.newZoom);
    }


    public void doIncrementZoom(final int zoomIncrement, final double newCenter) {
          doSetZoomCenter(getZoom() + zoomIncrement, newCenter);
    }

    /**
     * Intended to be called by UI elements, this method
     * performs all actions necessary to set a new zoom
     * and center location
     * @param newZoom
     * @param newCenter
     */
    public void doSetZoomCenter(final int newZoom, final double newCenter) {

        if (chrName.equals(Globals.CHR_ALL)) {
            chrName = getGenome().getHomeChromosome();
        }

        if (chrName.equals(Globals.CHR_ALL)) {
            // Translate the location to chromosome number
            synchronized (this){
                jumpToChromosomeForGenomeLocation(newCenter);
            }
            IGV.getInstance().chromosomeChangeEvent(chrName);
        } else {
            setZoom(newZoom);
            // Adjust origin so newCenter is centered
            centerOnLocation(newCenter);
        }

        IGV.repaintPanelsHeadlessSafe();
        recordHistory();
    }

    protected double getGenomeCenterPosition() {
        return origin + ((widthInPixels / 2) * getScale());
    }

    /**
     * Return the current locationScale in base pairs / pixel
     *
     * @return
     */
    public double getScale() {
        if (locationScale <= 0) {
            computeLocationScale();
        }
        return locationScale;
    }

    public void invalidateLocationScale() {
        this.locationScale = -1;
    }

    /**
     * Change the frame to the specified chromosome.
     *
     * @param name
     * @param force
     */
    public synchronized void setChromosomeName(String name, boolean force) {

        if ((chrName == null) || !name.equals(chrName) || force) {
            chrName = name;
            origin = 0;
            setEnd = -1;

            this.zoom = -1;
            setZoom(0);

            chromoObservable.setChangedAndNotify();
        }
    }

    /**
     * Recalculate the locationScale, based on {@link #setEnd}, {@link #origin}, and
     * {@link #widthInPixels}
     * DOES NOT alter zoom value
     */
    private synchronized void computeLocationScale() {
        Genome genome = getGenome();

        //Should consider getting rid of this. We don't have
        //a chromosome length without a genome, not always a problem
        if (genome != null) {
            if (setEnd > 0 && widthInPixels > 0) {
                this.locationScale = ((setEnd - origin) / widthInPixels);
                setEnd = -1;
            } else {
                double virtualPixelSize = getTilesTimesBinsPerTile();
                double nPixel = Math.max(virtualPixelSize, widthInPixels);
                this.locationScale = (((double) getChromosomeLength()) / nPixel);
            }
        }
    }

    /**
     * Recalculate the zoom value based on current start/end
     * locationScale is not altered
     *
     */
    private void computeZoom() {
        int newZoom = calculateZoom(getOrigin(), getEnd());
        setZoomWithinLimits(newZoom);
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
        }
    }

    /**
     * Record the current state of the frame in history.
     * It is recommended that this NOT be called from within ReferenceFrame,
     * and callers use it after making all changes
     */
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
        IGV.repaintPanelsHeadlessSafe();
    }

    public void snapToGrid() {
        setOrigin(Math.round(origin));
        IGV.repaintPanelsHeadlessSafe();
    }

    public void centerOnLocation(String chr, double chrLocation) {
        if (!chrName.equals(chr)) {
            setChromosomeName(chr);
        }
        centerOnLocation(chrLocation);
    }

    public void centerOnLocation(double chrLocation) {
        double windowWidth = (widthInPixels * getScale()) / 2;
        setOrigin(Math.round(chrLocation - windowWidth));
        IGV.repaintPanelsHeadlessSafe();
    }

    public boolean windowAtEnd() {
        double windowLengthBP = widthInPixels * getScale();
        return origin + windowLengthBP + 1 > getChromosomeLength();
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
            if (start >= 0 && end >= 0) {
                this.origin = start;
                this.setEnd = end;
                computeLocationScale();
                setZoomWithinLimits(calculateZoom(this.getOrigin(), this.getEnd()));
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


    /**
     * Compute the maximum zoom level, which is a function of chromosome length.
     */
    private void computeMaxZoom() {
        Genome genome = getGenome();
        // Compute max zoom.  Assume window size @ max zoom of ~ 50 bp
        if (genome != null && chrName != null && genome.getChromosome(chrName) != null) {
            if (chrName.equals(Globals.CHR_ALL)) {
                maxZoom = 0;
            } else {
                int chrLength = genome.getChromosome(chrName).getLength();
                maxZoom = (int) (Math.log(chrLength / 50.0) / Globals.log2) + 1;
            }
        }
    }

    /**
     * Calculate the zoom level given start/end in bp.
     * Doesn't change anything
     * @param start
     * @param end
     * @return
     */
    protected int calculateZoom(double start, double end) {
        return (int) (Math.log(getChromosomeLength() / (end - start)) / Globals.log2) + 1;
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
        return getTilesTimesBinsPerTile();
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
    public int getWidthInPixels() {
        return widthInPixels;
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
        return nTiles * (double) binsPerTile;
    }

    public int getMidpoint() {
        return pixelX + widthInPixels / 2;
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
