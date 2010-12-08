/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.Locus;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;

import java.text.NumberFormat;

/**
 * @author jrobinso
 */
public class ReferenceFrame {

    private static Logger log = Logger.getLogger(ReferenceFrame.class);

    String name;
    int pixelX;
    private int pixelWidth;

    /**
     * The nominal viewport width in pixels.
     */
    public static int binsPerTile = 700;

    /**
     * The chromosome currently in view
     */
    private String chrName = "chrAll";

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
    private double origin = 0;

    private double end = -1;

    /**
     * The location (x axis) locationScale in base pairs / virtual pixel
     */
    private double locationScale;

    private boolean locationScaleValid = false;


    public ReferenceFrame(String name) {
        this.name = name;
        Genome genome = getGenome();
        this.chrName = genome == null ? "" : genome.getHomeChromosome();
    }


    public void setBounds(int x, int w) {
        this.pixelX = x;
        if (w != pixelWidth) {
            pixelWidth = w;
            computeLocationScale();
        }
    }

    public int getMidpoint() {
        return pixelX + pixelWidth / 2;
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
                setZoom(maxZoom);
            }
        }


        // TODO -- fire chromosome changed event rather than this
    }

    private void setZoom(int newZoom) {

        // All  zoom events "release" the frame, enabling pan and zoom
        end = -1;   // <= A lousy convention, means end is computed (free to float)

        zoom = Math.min(maxZoom, newZoom);
        nTiles = (int) Math.pow(2, Math.max(minZoom, zoom));
        maxPixel = getTilesTimesBinsPerTile();
        invalidateLocationScale();

        // TODO -- do this with events,
        if (IGVMainFrame.hasInstance()) {
            IGVMainFrame.getInstance().repaintStatusAndZoomSlider();
        }

    }

    public void incrementZoom(int increment) {
        zoomAndCenter(zoom + increment);
    }

    public void zoomAndCenterAdjusted(int newZoom) {
        zoomAndCenter(minZoom + newZoom);
    }

    public void zoomAndCenter(int newZoom) {

        // Zoom but remain centered about current center
        double currentCenter = origin + ((pixelWidth / 2) * getScale());

        zoomTo(newZoom, currentCenter);
    }

    public void zoomTo(final int newZoom, final double newCenter) {

        if (chrName.equals(Globals.CHR_ALL)) {
            chrName = getGenome().getHomeChromosome();
        }

        if (chrName.equals(Globals.CHR_ALL)) {

            // DISABLE ZOOMING FOR GENOME VIEW
            // Translate the location to chromosome number
            jumpToChromosomeForGenomeLocation(newCenter);
            IGVMainFrame.getInstance().chromosomeChangeEvent();
        } else {
            if (zoom != newZoom) {

                setZoom(newZoom);
                computeLocationScale();

                double newLocationScale = getScale();

                // Adjust origin so newCenter is centered
                double newOrigin = Math.round(
                        newCenter - ((pixelWidth / 2) * newLocationScale));

                setOrigin(newOrigin);

            }
        }

        recordHistory();

        // This is a hack,  this is not a drag event but is a "jump"
        DragEventManager.getInstance().dragStopped();

    }

    public void recordHistory() {
        IGVMainFrame.getInstance().getSession().getHistory().push(getCurrentLocusString());
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
                setZoom(maxZoom);
            }
        }
        double windowWidth = (pixelWidth * getScale()) / 2;
        setOrigin(Math.round(chrLocation - windowWidth));
    }

    public void centerOnLocation(double chrLocation) {
        double windowWidth = (pixelWidth * getScale()) / 2;
        setOrigin(Math.round(chrLocation - windowWidth));
        recordHistory();
    }

    public boolean windowAtEnd() {
        double windowLengthBP = pixelWidth * getScale();
        return origin + windowLengthBP + 1 > getChromosomeLength();

    }


    /* Keep origin within data range */

    /**
     * Method description
     *
     * @param newOrigin
     */
    public void setOrigin(double newOrigin) {

        // All movements "release" the frame, enabling pan and zoom
        //end = -1;   // <= A lousy convention, means end is computed (free to float)

        int windowLengthBP = (int) (pixelWidth * getScale());
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_SOFT_CLIPPED)) {
            origin = Math.max(-1000, Math.min(newOrigin, getChromosomeLength() + 1000 - windowLengthBP));
        } else {
            origin = Math.max(0, Math.min(newOrigin, getChromosomeLength() - windowLengthBP));
        }

        // If zoomed in sufficiently track the center position
        if (locationScale < 10) {
            IGVMainFrame.getInstance().setStatusBarMessage(chrName + ":" + ((int) getCenter() + 1));
        }

        // Repaint
        IGVMainFrame.getInstance().repaintDataAndHeaderPanels();
        IGVMainFrame.getInstance().repaintStatusAndZoomSlider();
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

            if (pixelWidth > 0) {
                setLocationScale(((double) (end - start)) / pixelWidth);
            }
            origin = start;
            this.end = end;
        }

        if (log.isDebugEnabled()) {
            log.debug("Data panel width = " + pixelWidth);
            log.debug("New start = " + (int) origin);
            log.debug("New end = " + (int) getEnd());
            log.debug("New center = " + (int) getCenter());
            log.debug("Scale = " + locationScale);
        }


        // Repaint
        IGVMainFrame.getInstance().repaintDataAndHeaderPanels();
        IGVMainFrame.getInstance().repaintStatusAndZoomSlider();

    }

    private void imputeZoom(double start, double end) {
        int z = (int) (Math.log(getChromosomeLength() / (end - start)) / Globals.log2) + 1;
        if (z != this.zoom) {
            zoom = Math.min(maxZoom, Math.max(minZoom, z));
            nTiles = (int) Math.pow(2, zoom);
            maxPixel = getTilesTimesBinsPerTile();
        }
    }

    private Genome getGenome() {
        return GenomeManager.getInstance().getCurrentGenome();
    }

    public double getOrigin() {
        return origin;
    }

    public double getCenter() {
        return origin + getScale() * pixelWidth / 2;
    }

    public double getEnd() {
        return origin + getScale() * pixelWidth;
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
            setZoom(0);
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

    public int getPixelWidth() {
        return pixelWidth;
    }

    /**
     * Return the current locationScale in base pairs / pixel
     *
     * @return
     */
    public double getScale() {
        if ((locationScale == 0) || !locationScaleValid) {
            computeLocationScale();
        }

        if (locationScale < 0) {
            System.err.println("Negative scale");
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

            computeMinZoom();

            if (end < 0) {
                double virtualPixelSize = getTilesTimesBinsPerTile();

                double nPixel = Math.max(virtualPixelSize, pixelWidth);

                setLocationScale(((double) getChromosomeLength()) / nPixel);
            } else {
                int w = pixelWidth;
                setLocationScale((end - origin) / w);
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
            minZoom = Math.max(0, (int) (Math.log((pixelWidth / binsPerTile)) / Globals.log2));

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
            log.error("Genome not loaded!");
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
                if (genome == null) {
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
    public String getCurrentLocusString() {

        if (zoom == 0) {
            return getChrName();
        } else {

            Range range = getCurrentRange();
            String startStr = NumberFormat.getInstance().format(range.getStart());
            String endStr = NumberFormat.getInstance().format(range.getEnd());
            String position = range.getChr() + ":" + startStr + "-" + endStr;

            return position;
        }
    }

    public Range getCurrentRange() {
        int start = 0;
        int end = pixelWidth;
        int startLoc = (int) getChromosomePosition(start) + 1;
        int endLoc = (int) getChromosomePosition(end);
        Range range = new Range(getChrName(), startLoc, endLoc);
        return range;
    }

    private void setLocationScale(double locationScale) {
        System.out.println("Set location scale: " + locationScale + "  origin=" + origin + "  end=" + end);
        this.locationScale = locationScale;
        locationScaleValid = true;

    }

    public void setInterval(Locus locus) {
        setInterval(locus.getChr(), locus.getStart(), locus.getEnd());
    }


    public void setInterval(String chr, int start, int end) {
        this.chrName = chr;
        this.origin = start;
        this.end = end;
        if (pixelWidth > 0) {
            locationScale = (end - origin) / pixelWidth;
            locationScaleValid = true;
        }
        imputeZoom(origin, end);
    }


    public void reset() {
        setInterval(FrameManager.getLocus(name));
    }

    public String getName() {
        return name;
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
