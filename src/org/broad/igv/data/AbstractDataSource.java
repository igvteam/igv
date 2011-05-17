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
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.LRUCache;

import java.util.*;

/**
 * @author jrobinso
 */
public abstract class AbstractDataSource implements DataSource {

    private static Logger log = Logger.getLogger(AbstractDataSource.class);

    // DataManager dataManager;
    boolean cacheSummaryTiles = true;
    WindowFunction windowFunction = WindowFunction.mean;
    LRUCache<String, SummaryTile> summaryTileCache = new LRUCache(this, 20);
    Genome genome;

    public AbstractDataSource(Genome genome) {
        this.genome = genome;
    }

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
    protected List<SummaryTile> getPrecomputedSummaryTiles(String chr, int startLocation, int endLocation, int zoom) {
        return null;
    }


    public int getChrLength(String chr) {
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
    //        String genomeId = IGV.getInstance().getGenomeManager().getGenomeId();
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
        //FeatureUtils.sortFeatureList(summaryScores);
        return summaryScores;

    }

    private List<SummaryTile> getSummaryTilesForRange(String chr, int startLocation, int endLocation, int zReq) {
        int startTile;
        int endTile;
        if (zReq < getNumZoomLevels(chr)) {
            return getPrecomputedSummaryTiles(chr, startLocation, endLocation, zReq);
        } else {
            int chrLength = getChrLength(chr);
            if (chrLength == 0) {
                return Collections.emptyList();
            }
            endLocation = Math.min(endLocation, chrLength);

            // By definition there are 2^z tiles per chromosome, and 700 bins per tile, where z is the zoom level.
            //int maxZoom = (int) (Math.log(chrLength/700) / Globals.log2) + 1;
            //int z = Math.min(zReq, maxZoom);
            int z = zReq;
            int nTiles = (int) Math.pow(2, z);
            double binSize = Math.max(1, (((double) chrLength) / nTiles) / 700);



            int adjustedStart = Math.max(0, startLocation);
            int adjustedEnd = Math.min(chrLength, endLocation);


            if (cacheSummaryTiles && !FrameManager.isGeneListMode()) {
                double tileWidth = ((double) chrLength) / nTiles;

                startTile = (int) (adjustedStart / tileWidth);
                endTile = (int) (Math.min(chrLength, adjustedEnd) / tileWidth) + 1;
                List<SummaryTile> tiles = new ArrayList(nTiles);
                for (int t = startTile; t <= endTile; t++) {
                    int tileStart = (int) (t * tileWidth);
                    int tileEnd = Math.min(chrLength, (int) ((t + 1) * tileWidth));

                    String key = chr + "_" + z + "_" + t + getWindowFunction();
                    SummaryTile summaryTile = summaryTileCache.get(key);
                    if (summaryTile == null) {

                        summaryTile = computeSummaryTile(chr, t, tileStart, tileEnd, binSize);

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
            } else {
                SummaryTile summaryTile = computeSummaryTile(chr, 0, startLocation, endLocation, binSize);
                return Arrays.asList(summaryTile);
            }


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

    SummaryTile computeSummaryTile(String chr, int tileNumber, int startLocation, int endLocation, double binSize) {


        // TODO -- we should use an index here
        int longestGene = getLongestFeature(chr);

        int adjustedStart = Math.max(startLocation - longestGene, 0);
        DataTile rawTile = getRawData(chr, adjustedStart, endLocation);
        SummaryTile tile = null;

        int nBins = (int) ((endLocation - startLocation) / binSize + 1);

        if ((rawTile != null) && !rawTile.isEmpty() && nBins > 0) {
            int[] starts = rawTile.getStartLocations();
            int[] ends = rawTile.getEndLocations();
            float[] values = rawTile.getValues();
            String[] features = rawTile.getFeatureNames();

            tile = new SummaryTile(tileNumber, startLocation);

            if (windowFunction == WindowFunction.none || FrameManager.isGeneListMode()) {

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

                    BasicScore score = new BasicScore(s, e, v);
                    tile.addScore(score);

                }


            } else {
                Bin[] bins = new Bin[nBins];
                
                List<LocusScore> scores = new ArrayList(nBins);

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
                        if (b < bins.length) {
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

                }


                // Aggregate adjacent bins.  This stiches back together features that span multiple bins.
                // TODO-- look at computing variable length bins to start with

                /*Bin currentBin = null;
                for (int b = 0; b < bins.length; b++) {
                    if (bins[b] != null) {
                        if (b < bins.length) {
                            if (currentBin == null) {
                                currentBin = bins[b];
                            } else {
                                if (false) { ///currentBin.isExtension(bins[b])) {
                                    currentBin.setEnd(bins[b].getEnd());
                                } else {
                                    scores.add(currentBin);
                                    currentBin = bins[b];
                                }
                            }
                        }
                    }
                }
                if (currentBin != null) {
                    scores.add(currentBin);
                }*/
                for(Bin bin : bins) {
                    if(bin != null)
                    scores.add(bin);
                }

                tile.addAllScores(scores);
            }
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
        this.summaryTileCache.clear();
    }


    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    // TODO -- get window functions dynamically from data
    static List<WindowFunction> wfs = new ArrayList();


    static {
        wfs.add(WindowFunction.min);
        wfs.add(WindowFunction.percentile10);
        wfs.add(WindowFunction.median);
        wfs.add(WindowFunction.mean);
        wfs.add(WindowFunction.percentile90);
        wfs.add(WindowFunction.max);

    }


    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return wfs;
    }



    /*
    SummaryTile computeSummaryTile(String chr, int tileNumber, int startLocation, int endLocation, float binSize) {


        // TODO -- we should use an index here
        int longestGene = getLongestFeature(chr);

        int adjustedStart = Math.max(startLocation - longestGene, 0);
        DataTile rawTile = getRawData(chr, adjustedStart, endLocation);
        SummaryTile tile = null;

        int nBins = (int) ((endLocation - startLocation) / binSize + 1);

        if ((rawTile != null) && !rawTile.isEmpty() && nBins > 0) {
            int[] starts = rawTile.getStartLocations();
            int[] ends = rawTile.getEndLocations();
            float[] values = rawTile.getValues();
            String[] features = rawTile.getFeatureNames();

            tile = new SummaryTile(tileNumber, startLocation);

            if (windowFunction == WindowFunction.none) {

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

                    BasicScore score = new BasicScore(s, e, v);
                    tile.addScore(score);

                }


            } else {
                List<Bin> bins = new ArrayList();

                int startIdx = 0;
                int endIdx = starts.length;
                for (int i = 0; i < starts.length; i++) {
                    int s = starts[i];
                    int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);
                    if (e < startLocation) {
                        startIdx = i;
                        continue;
                    } else if (s > endLocation) {
                        endIdx = i;
                        break;
                    }
                }


                for (int i = startIdx; i < endIdx; i++) {
                    int s = starts[i];
                    int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);
                    float v = values[i];
                    String probeName = features == null ? null : features[i];

                    int start = binPosition(binSize, s);
                    int end = binPosition(binSize, e);

                    // Previous bins are already sliced into non-overlapping section.  The current bin can (1) contribute
                    // to previous bins, (2) slice a previous bin if end intersects it, or (3) start a new bin
                    // if start is >= the end of the last previous bin

                    int lastEnd = 0;
                    int binIdx = 0;
                    for (Bin b : bins) {
                        int bEnd = b.getEnd();
                        lastEnd = bEnd;
                        if (b.getStart() >= end) {
                            // This should never happen?
                            break;
                        }
                        if (b.getEnd() <= start) {
                            continue;
                        }

                        if (b.getStart() < start) {
                            //    [     ]            <= b
                            //        [     ]
                            // Shift end of bin
                            b.setEnd(start);

                            // Create the merged bin, representing the overlapped region
                            Bin mergedBin = new Bin(start, bEnd, probeName, v, windowFunction);
                            mergedBin.mergeValues(b);
                            bins.add(mergedBin);

                            // Create a new bin for the non-overlapped part
                            Bin newBin = new Bin(bEnd, end, probeName, v, windowFunction);
                            bins.add(newBin);

                        } else {
                            //     [     ]          <= b
                            //       [ ]
                            // Shift end of bin
                            b.setEnd(start);

                            // Create the merged bin, representing the overlapped region
                            Bin mergedBin = new Bin(start, bEnd, probeName, v, windowFunction);
                            mergedBin.mergeValues(b);
                            bins.add(mergedBin);

                            // Create a new bin for the non-overlapped part
                            Bin newBin = new Bin(bEnd, end, windowFunction);
                            newBin.mergeValues(b);
                            bins.add(newBin);
                        }
                        binIdx++;
                    }

                    if (start >= lastEnd) {
                        bins.add(new Bin(start, end, probeName, v, windowFunction));
                    }

                }

                tile.addAllScores(bins);
            }
        }

        return tile;
    }

    private int binPosition(float binSize, int s) {
        return (int) (Math.round(s / binSize) * binSize);
    }  */
    
}
