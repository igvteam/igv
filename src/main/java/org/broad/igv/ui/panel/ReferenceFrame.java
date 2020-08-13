/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.event.ShiftEvent;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.InsertionManager;
import org.broad.igv.sam.InsertionMarker;
import org.broad.igv.ui.IGV;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.ViewChange;
import org.broad.igv.ui.util.MessageUtils;


/**
 * @author jrobinso
 */
public class ReferenceFrame {

    private static Logger log = Logger.getLogger(ReferenceFrame.class);

    IGVEventBus eventBus;

    /**
     * The origin in bp
     */
    public volatile double origin = 0;

    /**
     * The nominal viewport width in pixels.
     */
    public static int binsPerTile = 700;

    boolean visible = true;

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
     * Minimum allowed range in base-pairs
     */
    protected static final int minBP = 40;

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
    protected double nTiles = 1;


    /**
     * The location (x axis) locationScale in base pairs / virtual pixel
     */
    protected volatile double scale;

    protected Locus initialLocus = null;


    public ReferenceFrame(String name) {
        this.name = name;
        Genome genome = getGenome();
        this.chrName = genome == null ? "" : genome.getHomeChromosome();
        this.eventBus = IGVEventBus.getInstance();
    }

    public ReferenceFrame(ReferenceFrame otherFrame) {
        this(otherFrame, otherFrame.eventBus);
    }

    /**
     * Copy constructor with event bus ovverride -- used by Sashimii plot
     *
     * @param otherFrame
     */
    public ReferenceFrame(ReferenceFrame otherFrame, IGVEventBus eventBus) {
        this.chrName = otherFrame.chrName;
        this.initialLocus = otherFrame.initialLocus;
        this.scale = otherFrame.scale;
        this.minZoom = otherFrame.minZoom;
        this.name = otherFrame.name;
        this.nTiles = otherFrame.nTiles;
        this.origin = otherFrame.origin;
        this.pixelX = otherFrame.pixelX;
        this.widthInPixels = otherFrame.widthInPixels;
        this.zoom = otherFrame.zoom;
        this.maxZoom = otherFrame.maxZoom;
        this.eventBus = eventBus;
    }

    public boolean isVisible() {
        return visible;
    }

    public void setVisible(boolean visible) {
        this.visible = visible;
    }


    public void dragStopped() {
        setOrigin(Math.round(origin));   // Snap to gride
        eventBus.post(ViewChange.Result());
    }

    public void changeGenome(Genome genome) {
        setChromosomeName(genome.getHomeChromosome(), true);
    }

    public void changeChromosome(String chrName, boolean recordHistory) {
        boolean changed = setChromosomeName(chrName, false);
        // if (changed) {
        ViewChange resultEvent = ViewChange.ChromosomeChangeResult(chrName);
        resultEvent.setRecordHistory(recordHistory);
        eventBus.post(resultEvent);
        changeZoom(0);
        // }
    }

    public void changeZoom(int newZoom) {
        doSetZoom(newZoom);
        ViewChange result = ViewChange.Result();
        result.setRecordHistory(false);
        eventBus.post(result);
    }


    /**
     * Set the position and width of the frame, in pixels
     * The origin/end positions are kept fixed iff valid
     *
     * @param pixelX
     * @param widthInPixels
     */
    public void setBounds(int pixelX, int widthInPixels) {
        this.pixelX = pixelX;

        if (this.widthInPixels != widthInPixels) {

            //If we have what looks like a valid end position we keep it
            if (this.widthInPixels > 0 && this.initialLocus == null) {
                int start = (int) getOrigin();
                int end = (int) getEnd();
                if (start >= 0 && end >= 1) {
                    this.initialLocus = new Locus(getChrName(), start, end);
                }
            }

            this.widthInPixels = widthInPixels;
            computeLocationScale();
            computeZoom();
        }

    }

    /**
     * Sets zoom level and recomputes scale, iff newZoom != oldZoom
     * min/maxZoom are recalculated and respected,
     * and the locationScale is recomputed
     *
     * @param newZoom
     */
    protected void setZoom(int newZoom) {
        if (zoom != newZoom) {
            synchronized (this) {
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
        final IGVPreferences preferences = PreferencesManager.getPreferences();
        if (preferences.getAsBoolean(Constants.SAM_SHOW_SOFT_CLIPPED)) {
            newOrigin = Math.max(
                    -preferences.getAsInt(Constants.SAM_MAX_SOFT_CLIP),
                    Math.min(position, getMaxCoordinate() + preferences.getAsInt(Constants.SAM_MAX_SOFT_CLIP) - windowLengthBP));
        } else {
            newOrigin = Math.max(0, Math.min(position, getMaxCoordinate() - windowLengthBP));
        }
        origin = newOrigin;
    }


    protected synchronized void setZoomWithinLimits(int newZoom) {
        zoom = Math.max(minZoom, Math.min(maxZoom, newZoom));
        nTiles = Math.pow(2, zoom);
    }

    /**
     * Increment the zoom level by {@code zoomIncrement}, leaving
     * the center the same
     *
     * @param zoomIncrement
     */
    public void doZoomIncrement(int zoomIncrement) {
        double currentCenter = getGenomeCenterPosition();
        doIncrementZoom(zoomIncrement, currentCenter);
    }

    /**
     * Set the zoom level to {@code newZoom}, leaving
     * the center the same
     *
     * @param newZoom
     */
    public void doSetZoom(int newZoom) {
        double currentCenter = getGenomeCenterPosition();
        doSetZoomCenter(newZoom, currentCenter);
    }


    public void doIncrementZoom(final int zoomIncrement, final double newCenter) {
        doSetZoomCenter(getZoom() + zoomIncrement, newCenter);
    }

    /**
     * Intended to be called by UI elements, this method
     * performs all actions necessary to set a new zoom
     * and center location
     *
     * @param newZoom
     * @param newCenter Center position, in genome coordinates
     */
    public void doSetZoomCenter(final int newZoom, final double newCenter) {

        if (chrName.equals(Globals.CHR_ALL)) {
            chrName = getGenome().getHomeChromosome();
        }

        if (!chrName.equals(Globals.CHR_ALL)) {
            setZoom(newZoom);
            // Adjust origin so newCenter is centered
            centerOnLocation(newCenter);
        }
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
        if (scale <= 0) {
            computeLocationScale();
        }
        return scale;
    }

    /**
     * Calls {@link #setChromosomeName(String, boolean)} with force = false
     * It is preferred that you post an event to the EventBus instead, this is public
     * as an implementation side effect
     *
     * @param name
     * @return boolean indicating whether the chromosome actually changed
     */
    public boolean setChromosomeName(String name) {
        return setChromosomeName(name, false);
    }

    /**
     * Change the frame to the specified chromosome, clearing all
     * view parameters (zoom, locationScale) in the process
     *
     * @param name  Name of the new chromosome
     * @param force Whether to force a change to the new chromosome, even if it's
     *              the same name as the old one
     * @return boolean indicating whether the chromosome actually changed
     */
    public synchronized boolean setChromosomeName(String name, boolean force) {

        if (shouldChangeChromosome(name) || force) {
            chrName = name;
            origin = 0;
            this.scale = -1;
            this.calculateMaxZoom();

            this.zoom = -1;
            setZoom(0);

            //chromoObservable.setChangedAndNotify();
            return true;
        }

        return false;
    }


    /**
     * Record the current state of the frame in history.
     * It is recommended that this NOT be called from within ReferenceFrame,
     * and callers use it after making all changes
     * <p>
     * //TODO Should we save history by receiving events in History?
     */
    public void recordHistory() {
        IGV.getInstance().getSession().getHistory().push(getFormattedLocusString(), zoom);
    }


    public void shiftOriginPixels(int delta) {

        double shiftBP = delta * getScale();
        setOrigin(origin + shiftBP);
        eventBus.post(ViewChange.Result());
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
        eventBus.post(ViewChange.LocusChangeResult(chrName, origin, chrLocation + windowWidth));
    }

    public boolean windowAtEnd() {
        double windowLengthBP = widthInPixels * getScale();
        return origin + windowLengthBP + 1 > getMaxCoordinate();
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

        end = Math.min(getMaxCoordinate(chr), end);

 //       synchronized (this) {
            this.initialLocus = locus;
            this.chrName = chr;
            if (start >= 0 && end >= 0) {
                this.origin = start;
                beforeScaleZoom(locus);
                computeLocationScale();
                computeZoom();
            }
   //     }

        if (log.isDebugEnabled()) {
            log.debug("Data panel width = " + widthInPixels);
            log.debug("New start = " + (int) origin);
            log.debug("New end = " + (int) getEnd());
            log.debug("New center = " + (int) getCenter());
            log.debug("Scale = " + scale);
        }

        eventBus.post(ViewChange.LocusChangeResult(chrName, start, end));
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

    public void setAdjustedZoom(int zoom) {
        this.doSetZoom(minZoom + zoom);
    }

    /**
     * Determine if this view will change at all based on the {@code newChrName}
     * The view changes if newChrName != {@code #this.chr} or if we are not
     * at full chromosome view
     *
     * @param newChrName
     * @return
     */
    private boolean shouldChangeChromosome(String newChrName) {
        return chrName == null || !chrName.equals(newChrName);
    }


    protected void calculateMaxZoom() {
        this.maxZoom = Globals.CHR_ALL.equals(this.chrName) ? 0 :
                (int) Math.ceil(Math.log(getChromosomeLength() / minBP) / Globals.log2);
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

        InsertionMarker i = InsertionManager.getInstance().getSelectedInsertion(getChrName());

        if (i != null && i.position > origin) {
            // if (IGV.getInstance().getSession().expandInsertions && insertionMarkers != null && insertionMarkers.size() > 0) {
            double start = getOrigin();
            double scale = getScale();
            double iEnd = 0,
                    iStart = 0;


            iStart = iEnd + (i.position - start) / scale; // Screen position of insertionMarker start
            if (screenPosition < iStart) {
                return start + scale * (screenPosition - iEnd);
            }

            iEnd = iStart + i.size / scale;  // Screen position of insertionMarker end
            if (screenPosition < iEnd) {
                return i.position;   // In the gap
            }

            start = i.position + 1;
            //    }
            return start + scale * (screenPosition - iEnd);

        } else {
            return origin + getScale() * screenPosition;
        }

    }


    /**
     * Return the screen position corresponding to the chromosomal position.
     *
     * @param chromosomePosition
     * @return
     */
    public int getScreenPosition(double chromosomePosition) {

        InsertionMarker i = InsertionManager.getInstance().getSelectedInsertion(chrName);

        if (i == null || i.position < origin || i.position > chromosomePosition) {
            return (int) ((chromosomePosition - origin) / getScale());
        } else {
            return (int) ((chromosomePosition + i.size - origin) / getScale());
        }
    }


    public Chromosome getChromosome() {
        Genome genome = getGenome();
        if (genome == null) {
            return null;
        }

        return genome.getChromosome(chrName);
    }

    /**
     * The maximum coordinate currently allowed.
     * In genomic coordinates this is the same as the chromosome length.
     * In exome coordinates, the two are different
     * (since ExomeReferenceFrame takes input in genomic coordinates)
     *
     * @return
     * @see #getChromosomeLength()
     */
    public int getMaxCoordinate() {
        return this.getChromosomeLength();
    }

    private static int getMaxCoordinate(String chrName) {
        return getChromosomeLength(chrName);
    }

    /**
     * Chromosome length, in genomic coordinates.
     * Intended to be used for scaling
     *
     * @return
     * @see #getMaxCoordinate()
     */
    public int getChromosomeLength() {
        return getChromosomeLength(this.chrName);
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
        int endLoc = (int) Math.round(getChromosomePosition(widthInPixels));
        Range range = new Range(getChrName(), (int) origin, endLoc);
        return range;
    }

    public void reset() {
        jumpTo(FrameManager.getLocus(name));
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    /**
     * Recalculate the locationScale, based on {@link #initialLocus}, {@link #origin}, and
     * {@link #widthInPixels}
     * DOES NOT alter zoom value
     */
    private synchronized void computeLocationScale() {
        Genome genome = getGenome();

        //Should consider getting rid of this. We don't have
        //a chromosome length without a genome, not always a problem
        if (genome != null) {

            // The end location, in base pairs.
            // If negative, we use the whole chromosome
            int setEnd = -1;
            if (this.initialLocus != null) setEnd = this.initialLocus.getEnd();

            if (setEnd > 0 && widthInPixels > 0) {
                this.scale = ((setEnd - origin) / widthInPixels);
                this.initialLocus = null;
            } else {
                double virtualPixelSize = getTilesTimesBinsPerTile();
                double nPixel = Math.max(virtualPixelSize, widthInPixels);
                this.scale = (((double) getChromosomeLength()) / nPixel);
            }
        }
    }

    /**
     * Recalculate the zoom value based on current start/end
     * locationScale is not altered
     */
    private void computeZoom() {
        int newZoom = calculateZoom(getOrigin(), getEnd());
        setZoomWithinLimits(newZoom);
    }


    /**
     * Called before scaling and zooming, during jumpTo.
     * Intended to be overridden
     *
     * @param locus
     */
    private void beforeScaleZoom(Locus locus) {
        calculateMaxZoom();
    }

    /**
     * Calculate the zoom level given start/end in bp.
     * Doesn't change anything
     *
     * @param start
     * @param end
     * @return
     */
    private int calculateZoom(double start, double end) {
        return (int) Math.round((Math.log((getChromosomeLength() / (end - start)) * (((double) widthInPixels) / binsPerTile)) / Globals.log2));
    }


    private static int getChromosomeLength(String chrName) {
        Genome genome = getGenome();

        if (genome == null) {
            return 1;
        }

        if (chrName.equals("All")) {
            // Genome coordinates are in kb => divde by 1000
            return (int) (genome.getNominalLength() / 1000);
        } else {
            Chromosome chromosome = genome.getChromosome(chrName);
            if (chromosome == null) {
                log.error("Null chromosome: " + chrName);
                if (genome.getChromosomes().size() == 0) {
                    return 1;
                } else {
                    return genome.getChromosomes().iterator().next().getLength();
                }
            }
            return chromosome.getLength();
        }
    }


    public int stateHash() {
        int result;
        long temp;
        temp = Double.doubleToLongBits(origin);
        result = (int) (temp ^ (temp >>> 32));
        result = 31 * result + chrName.hashCode();
        result = 31 * result + zoom;
        result = 31 * result + pixelX;
        result = 31 * result + widthInPixels;
        temp = Double.doubleToLongBits(scale);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }

    private static Genome getGenome() {
        return GenomeManager.getInstance().getCurrentGenome();
    }


}

