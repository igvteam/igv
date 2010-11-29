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
package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.tdf.Bin;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.LRUCache;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Collections;

/**
 * @author jrobinso
 */
public abstract class AbstractDataSource implements DataSource {

    private static Logger log = Logger.getLogger(AbstractDataSource.class);

    // DataManager dataManager;
    boolean cacheSummaryTiles = true;
    WindowFunction windowFunction = WindowFunction.mean;
    LRUCache<String, SummaryTile> summaryTileCache = new LRUCache(20);

    // abstract DataManager getDataManager();

    /**
     * Return the number of precomputed zoom levels for the given chromosome.
     *
     * @param chr
     * @return
     */
    abstract protected int getNumZoomLevels(String chr);

    // abstract protected TrackType getTrackType();

    /**
     * Return "raw" (i.e. not summarized) data for the specified interval.
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @return
     */
    abstract protected DataTile getRawData(String chr, int startLocation, int endLocation);

    /**
     * Return the precomputed summary tiles for the given locus and zoom level.  If
     * there are non return null.
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @param zoom
     * @return
     */
    protected List<SummaryTile> getPrecomputedSummaryTiles(String chr, int startLocation, int endLocation,
                                                           int zoom) {
        return null;
    }


    public int getChrLength(String chr) {
        final Genome genome = GenomeManager.getInstance().getGenome();
        if (chr.equals(Globals.CHR_ALL)) {
            return (int) (genome.getLength() / 1000);
        } else {
            Chromosome c = genome.getChromosome(chr);
            return c == null ? 0 : c.getLength();
        }
    }

    /**
     * Refresh the underlying data. Default implementation does nothing, subclasses
     * can override
     *
     * @param timestamp
     */
    public void refreshData(long timestamp) {

        // ignore --
    }

    /**
     * Return the longest feature in the dataset for the given chromosome.  This
     * is needed when computing summary data for a region.
     * <p/>
     * TODO - This default implementaiton is crude and should be overriden by subclasses.
     *
     * @param chr
     * @return
     */
    public abstract int getLongestFeature(String chr);

    //{
    //
    //    if (getTrackType() == TrackType.GENE_EXPRESSION) {
    //        String genomeId = GenomeManager.getInstance().getGenomeId();
    //        GeneManager gm = GeneManager.getGeneManager(genomeId);
    //        return (gm == null) ? 1000000 : gm.getLongestGeneLength(chr);
    //    } else {
    //        return 1000;
    //    }
    //}


    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {


        List<SummaryTile> tiles = getSummaryTilesForRange(chr, startLocation, endLocation, zoom);


        List<LocusScore> summaryScores = new ArrayList(tiles.size() * 700);

        for (SummaryTile tile : tiles) {
            summaryScores.addAll(tile.getScores());

        }
        FeatureUtils.sortFeatureList(summaryScores);
        return summaryScores;
    }

    private List<SummaryTile> getSummaryTilesForRange(String chr, int startLocation, int endLocation, int zoom) {
        assert endLocation >= startLocation;

        int startTile;
        int endTile;
        if (zoom < getNumZoomLevels(chr)) {
            return getPrecomputedSummaryTiles(chr, startLocation, endLocation, zoom);
        } else {
            int chrLength = getChrLength(chr);
            if (chrLength == 0) {
                return Collections.emptyList();
            }
            endLocation = Math.min(endLocation, chrLength);

            int adjustedStart = Math.max(0, startLocation);
            int adjustedEnd = Math.min(chrLength, endLocation);

            int z = zoom; //Math.min(8, zoom);
            int nTiles = (int) Math.pow(2, z);
            double tileWidth = ((double) chrLength) / nTiles;


            startTile = (int) (adjustedStart / tileWidth);
            endTile = (int) (Math.min(chrLength, adjustedEnd) / tileWidth) + 1;
            List<SummaryTile> tiles = new ArrayList(nTiles);
            for (int t = startTile; t <= endTile; t++) {
                int tileStart = (int) (t * tileWidth);
                int tileEnd = Math.min(chrLength, (int) ((t + 1) * tileWidth));

                String key = chr + "_" + zoom + "_" + t + getWindowFunction();
                SummaryTile summaryTile = summaryTileCache.get(key);
                if (summaryTile == null) {
                    summaryTile = computeSummaryTile(chr, t, tileStart, tileEnd);

                    if (cacheSummaryTiles && !FrameManager.isGeneListMode()) {
                        synchronized (summaryTileCache) {
                            summaryTileCache.put(key, summaryTile);
                        }
                    }
                }


                if (summaryTile != null) {
                    tiles.add(summaryTile);
                }
            }

            return tiles;
        }
    }

    /**
     * Note:  Package scope used so this method can be unit tested
     *
     * @param chr
     * @param tileNumber
     * @param startLocation
     * @param endLocation
     * @return
     */

    SummaryTile computeSummaryTile(String chr, int tileNumber, int startLocation, int endLocation) {


        // TODO -- we should use an index here
        int longestGene = getLongestFeature(chr);

        int adjustedStart = Math.max(startLocation - longestGene, 0);
        DataTile rawTile = getRawData(chr, adjustedStart, endLocation);
        SummaryTile tile = null;

        int nBins = Math.min(700, endLocation - startLocation);
        if ((rawTile != null) && !rawTile.isEmpty() && nBins > 0) {

            tile = new SummaryTile(tileNumber, startLocation);
            double binSize = ((float) (endLocation - startLocation)) / nBins;
            Bin[] bins = new Bin[nBins];


            int[] starts = rawTile.getStartLocations();
            int[] ends = rawTile.getEndLocations();
            float[] values = rawTile.getValues();
            String[] features = rawTile.getFeatureNames();

            List<LocusScore> scores = new ArrayList(700);

            /* for (int i = 0; i < starts.length; i++) {
               int s = starts[i];
                int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);

                if (e < startLocation) {
                    continue;
                } else if (s > endLocation) {
                    break;
                }
               String probeName = features == null ? null : features[i];
               scores.add(new Bin(s, e, probeName, values[i], windowFunction));

           }
            */

            for (int i = 0; i < starts.length; i++) {
                int s = starts[i];
                int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);

                if (e < startLocation) {
                    continue;
                } else if (s > endLocation) {
                    break;
                }

                String probeName = features == null ? null : features[i];
                float v = values[i];


                int startBin = (int) Math.max(0, ((s - startLocation) / binSize));
                int endBin = (int) Math.min(bins.length, ((e - startLocation) / binSize));
                if (startBin == endBin) {
                    endBin++;
                }
                for (int b = startBin; b < endBin; b++) {
                    Bin bin = bins[b];
                    if (bin == null) {
                        int start = (int) (startLocation + b * binSize);
                        int end = (int) (startLocation + (b + 1) * binSize);
                        bins[b] = new Bin(start, end, probeName, v, windowFunction);
                    } else {
                        bin.addValue(probeName, v);
                    }
                }

            }


            // Aggregate adjacent bins.  This stiches back together features that span multiple bins.
            // TODO-- look at computing variable length bins to start with

            Bin currentBin = null;
            for (int b = 0; b < bins.length; b++) {
                if (bins[b] != null) {
                    if (currentBin == null) {
                        currentBin = bins[b];
                    } else {
                        if (currentBin.isExtension(bins[b])) {
                            currentBin.setEnd(bins[b].getEnd());
                        } else {
                            scores.add(currentBin);
                            currentBin = bins[b];
                        }
                    }
                }
            }
            if (currentBin != null) {
                scores.add(currentBin);
            }


            tile.addAllScores(scores);
        }

        return tile;
    }


    /**
     * Return true if the data has been log normalized.
     *
     * @return
     */
    public boolean isLogNormalized() {
        return true;
    }


    public void setWindowFunction(WindowFunction statType) {
        this.windowFunction = statType;
    }


    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    // TODO -- get window functions dynamically from data
    static List<WindowFunction> wfs = new ArrayList();


    static {
        wfs.add(WindowFunction.percentile10);
        wfs.add(WindowFunction.median);
        wfs.add(WindowFunction.mean);
        wfs.add(WindowFunction.percentile90);
        wfs.add(WindowFunction.max);

    }


    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return wfs;
    }
}
