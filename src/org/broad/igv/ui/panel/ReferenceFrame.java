/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.panel;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.xome.Block;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;

import java.text.NumberFormat;
import java.util.List;
import java.util.Locale;


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
    private int zoom = minZoom;
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


    public ReferenceFrame(String name) {
        this.name = name;
        Genome genome = getGenome();
        this.chrName = genome == null ? "" : genome.getHomeChromosome();
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
    }


    public void setBounds(int x, int w) {
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


        // TODO -- fire chromosome changed event rather than this
    }

    private void zoomTo(int newZoom) {
        // All  zoom events "release" the frame, enabling pan and zoom

        zoom = Math.min(maxZoom, newZoom);
        nTiles = (int) Math.pow(2, Math.max(minZoom, zoom));
        maxPixel = getTilesTimesBinsPerTile();
        invalidateLocationScale();

        // TODO -- do this with events,
        if (IGV.hasInstance()) {
            IGV.getInstance().repaintStatusAndZoomSlider();
        }
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
        zoomByAndCenter(increment);
    }

    public void zoomAndCenterAdjusted(int newZoom) {
        double currentCenter = origin + ((widthInPixels / 2) * getScale());
        zoomTo(newZoom, currentCenter);
    }

    public void zoomByAndCenter(int zoomIncrement) {
        // Zoom but remain centered about current center
        double currentCenter = origin + ((widthInPixels / 2) * getScale());
        zoomBy(zoomIncrement, currentCenter);
    }


    public void zoomBy(final int zoomFactor, final double newCenter) {

        if (FrameManager.isGeneListMode()) {
            double f = Math.pow(2.0, zoomFactor);
            locationScale = Math.max(minScale, locationScale / f);
            double newOrigin = Math.round(newCenter - ((widthInPixels / 2) * locationScale));
            setOrigin(newOrigin);
            imputeZoom(origin, setEnd);
        } else {
            int newZoom = Math.max(0, zoom + zoomFactor);
            zoomTo(newZoom, newCenter);
        }
    }

    public void zoomTo(final int newZoom, final double newCenter) {

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
                computeLocationScale();
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

        for (String chr : getGenome().getChromosomeNames()) {
            double endMB = startMB + getGenome().getChromosome(chr).getLength() / 1000.0;

            if ((locationMB > startMB) && (locationMB <= endMB)) {

                // this.jumpTo(chr, -1, -1);
                this.setChromosomeName(chr);
                break;
            }

            startMB = endMB;
        }
    }

    public void shiftOriginPixels(double delta) {

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
        double windowWidth = (widthInPixels * getScale()) / 2;
        setOrigin(Math.round(chrLocation - windowWidth));
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

    /* Keep origin within data range */
    public void setOrigin(double position, boolean repaint) {


        int windowLengthBP = (int) (widthInPixels * getScale());
        double newOrigin = origin;
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_SOFT_CLIPPED)) {
            newOrigin = Math.max(-1000, Math.min(position, getChromosomeLength() + 1000 - windowLengthBP));
        } else {
            newOrigin = Math.max(0, Math.min(position, getChromosomeLength() - windowLengthBP));
        }
        double delta = newOrigin - origin;
        origin = newOrigin;

        // If zoomed in sufficiently track the center position
        //if (locationScale < 10) {
        //    IGV.getInstance().setStatusBarMessage(chrName + ":" + ((int) getCenter() + 1));
        //}

        // Repaint
        if (repaint) {
            IGV.getInstance().repaintDataAndHeaderPanels();
            IGV.getInstance().repaintStatusAndZoomSlider();
        }
    }

    public void jumpTo(String chr, int start, int end) {
        Genome genome = getGenome();
        if (chr != null) {
            if (genome.getChromosome(chr) == null && !chr.contains(Globals.CHR_ALL)) {
                MessageUtils.showMessage(chr + " is not a valid chromosome.");
                return;
            }
            setChromosomeName(chr);
        }

        Chromosome chromosome = genome == null ? null : genome.getChromosome(chr);
        if (chromosome != null) {
            end = Math.min(chromosome.getLength(), end);
        }

        if (start >= 0) {

            imputeZoom(start, end);

            if (widthInPixels > 0) {
                setLocationScale(((double) (end - start)) / widthInPixels);
            }
            origin = start;
        }

        if (log.isDebugEnabled()) {
            log.debug("Data panel width = " + widthInPixels);
            log.debug("New start = " + (int) origin);
            log.debug("New end = " + (int) getEnd());
            log.debug("New center = " + (int) getCenter());
            log.debug("Scale = " + locationScale);
        }


        // Repaint
        IGV.getInstance().repaintDataAndHeaderPanels();
        IGV.getInstance().repaintStatusAndZoomSlider();

    }

    protected void imputeZoom(double start, double end) {
        int z = (int) (Math.log(getChromosomeLength() / (end - start)) / Globals.log2) + 1;
        if (z != this.zoom) {
            zoom = Math.min(maxZoom, Math.max(minZoom, z));
            nTiles = (int) Math.pow(2, zoom);
            maxPixel = getTilesTimesBinsPerTile();
        }
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


    public void setChrName(String name) {
        this.setChromosomeName(name);
    }

    /**
     * @param name
     * @param force
     * @ deprecated, replace with calls to setChrName();
     */
    public void setChromosomeName(String name, boolean force) {

        if ((chrName == null) || !name.equals(chrName) || force) {
            chrName = name;
            origin = 0;
            zoomTo(0);
            computeMaxZoom();

            //dhmay: if there's a RegionNavigatorDialog around, need to update it.
            //Todo: Properly, there should be a way to register listeners to chromosome changes.
            if (RegionNavigatorDialog.getActiveInstance() != null)
                RegionNavigatorDialog.getActiveInstance().updateChromosomeDisplayed();
        }
    }

    /**
     * @param name
     * @ deprecated, replace with calls to setChrName();
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
    public double getScale() {
        if ((locationScale == 0) || !locationScaleValid || locationScale < 0) {
            computeLocationScale();
        }

        return locationScale;
    }


    public void invalidateLocationScale() {
        //Thread.dumpStack();
        locationScaleValid = false;
    }

    public void computeLocationScale() {
        Genome genome = getGenome();

        if (genome != null) {

            if (setEnd > 0) {
                locationScale = (setEnd - origin) / widthInPixels;
                locationScaleValid = true;
                imputeZoom(origin, setEnd);
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
    private void computeMinZoom() {
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
            return (int) (genome.getLength() / 1000);

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

    private void setLocationScale(double locationScale) {
        if (log.isDebugEnabled()) {
            log.debug("Set location scale: " + locationScale + "  origin=" + origin);
        }
        this.locationScale = locationScale;
        locationScaleValid = true;

    }

    public void setInterval(Locus locus) {
        this.initialLocus = locus;

        this.chrName = locus.getChr();
        this.origin = locus.getStart();
        if (widthInPixels > 0) {
            locationScale = (locus.getEnd() - origin) / widthInPixels;
            locationScaleValid = true;
            imputeZoom(origin, locus.getEnd());
        } else {
            // Set end temporarily until scale can be calculated
            this.setEnd = locus.getEnd();
        }
    }


    public void reset() {
        setInterval(FrameManager.getLocus(name));
        IGV.getInstance().repaintDataAndHeaderPanels();
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
