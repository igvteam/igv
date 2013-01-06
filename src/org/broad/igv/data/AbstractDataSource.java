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
package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tdf.Accumulator;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.collections.LRUCache;

import java.util.*;

/**
 * @author jrobinso
 */
public abstract class AbstractDataSource implements DataSource {

    private static Logger log = Logger.getLogger(AbstractDataSource.class);

    // DataManager dataManager;
    boolean cacheSummaryTiles = true;
    WindowFunction windowFunction = WindowFunction.mean;
    LRUCache<String, SummaryTile> summaryTileCache = new LRUCache(this, 10);
    protected Genome genome;

    public AbstractDataSource(Genome genome) {
        this.genome = genome;
    }


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
    abstract protected List<LocusScore> getPrecomputedSummaryScores(String chr, int startLocation, int endLocation, int zoom);


    public int getChrLength(String chr) {
        if (chr.equals(Globals.CHR_ALL)) {
            return (int) (genome.getNominalLength() / 1000);
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


        List<LocusScore> scores = getPrecomputedSummaryScores(chr, startLocation, endLocation, zoom);
        if (scores != null) {
            return scores;
        }

        List<SummaryTile> tiles = getSummaryTilesForRange(chr, startLocation, endLocation, zoom);

        scores = new ArrayList(tiles.size() * 700);

        for (SummaryTile tile : tiles) {
            scores.addAll(tile.getScores());

        }
        //FeatureUtils.sortFeatureList(summaryScores);
        return scores;

    }

    private List<SummaryTile> getSummaryTilesForRange(String chr, int startLocation, int endLocation, int zReq) {

        int chrLength = getChrLength(chr);
        if (chrLength == 0) {
            return Collections.emptyList();
        }
        endLocation = Math.min(endLocation, chrLength);


        int adjustedStart = Math.max(0, startLocation);
        int adjustedEnd = Math.min(chrLength, endLocation);


        if (cacheSummaryTiles && !FrameManager.isGeneListMode() && !FrameManager.isExomeMode()) {

            // By definition there are 2^z tiles per chromosome, and 700 bins per tile, where z is the zoom level.
            //int maxZoom = (int) (Math.log(chrLength/700) / Globals.log2) + 1;
            //int z = Math.min(zReq, maxZoom);
            int z = zReq;
            int virtualTileCount = (int) Math.pow(2, z);

            double tileWidth = ((double) chrLength) / virtualTileCount;
            int startTile = (int) (adjustedStart / tileWidth);
            int endTile = (int) (Math.min(chrLength, adjustedEnd) / tileWidth) + 1;
            List<SummaryTile> tiles = null;

            tiles = new ArrayList(endTile - startTile + 1);
            for (int t = startTile; t <= endTile; t++) {
                int tileStart = (int) (t * tileWidth);
                int tileEnd = Math.min(chrLength, (int) ((t + 1) * tileWidth));

                String key = chr + "_" + z + "_" + t + getWindowFunction();
                SummaryTile summaryTile = summaryTileCache.get(key);
                if (summaryTile == null) {

                    summaryTile = computeSummaryTile(chr, t, tileStart, tileEnd, 700);

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
            SummaryTile summaryTile = computeSummaryTile(chr, 0, startLocation, endLocation, 700);
            return Arrays.asList(summaryTile);
        }


    }


    /**
     * Note:  Package scope used so this method can be unit tested
     *
     * @param chr
     * @param tileNumber
     * @param startLocation
     * @param endLocation
     * @param nBins
     * @return
     */

    SummaryTile computeSummaryTile(String chr, int tileNumber, int startLocation, int endLocation, int nBins) {


        // TODO -- we should use an index here
        int longestGene = getLongestFeature(chr);

        int adjustedStart = Math.max(startLocation - longestGene, 0);
        DataTile rawTile = getRawData(chr, adjustedStart, endLocation);
        SummaryTile tile = null;


        if ((rawTile != null) && !rawTile.isEmpty() && nBins > 0) {
            int[] starts = rawTile.getStartLocations();
            int[] ends = rawTile.getEndLocations();
            float[] values = rawTile.getValues();
            String[] features = rawTile.getFeatureNames();

            tile = new SummaryTile();

            if (windowFunction == WindowFunction.none) {

                for (int i = 0; i < starts.length; i++) {
                    int s = starts[i];
                    int e = ends == null ? s + 1 : Math.max(s + 1, ends[i]);

                    if (e < startLocation) {
                        continue;
                    } else if (s >= endLocation) {
                        break;
                    }

                    String probeName = features == null ? null : features[i];
                    float v = values[i];

                    BasicScore score = new NamedScore(s, e, v, probeName);
                    tile.addScore(score);

                }


            } else {
                float normalizationFactor = 1.0f;
                List<LocusScore> scores = new ArrayList(nBins);
                double scale = (double) (endLocation - startLocation) / nBins;

                Accumulator accumulator = new Accumulator(windowFunction, 5);
                int accumulatedStart = -1;
                int accumulatedEnd = -1;
                int lastEndBin = 0;

                int size = starts.length;

                // Loop through and bin scores for this interval.
                for (int i = 0; i < size; i++) {

                    int true_end = ends == null ? starts[i] + 1 : ends[i];

                    float v = values[i] * normalizationFactor;
                    if (starts[i] >= endLocation) {
                        break;  // We're beyond the end of the requested interval
                    } else if (true_end <= startLocation || Float.isNaN(v)) {
                        //Not yet to interval, or not a number
                        continue;
                    }

                    // Bound feature at interval, other "piece" will be in another tile.
                    int s = Math.max(startLocation, starts[i]);
                    int e = Math.min(endLocation, true_end);

                    String probeName = features == null ? null : features[i];

                    // Compute bin numbers, relative to start of this tile
                    int endBin = (int) ((e - startLocation) / scale);
                    int startBin = (int) ((s - startLocation) / scale);

                    // If this feature spans multiple bins, or extends beyond last end bin, record
                    if (endBin > lastEndBin || endBin > startBin) {
                        if (accumulator.hasData()) {
                            scores.add(getCompositeScore(accumulator, accumulatedStart, accumulatedEnd));
                            accumulator = new Accumulator(windowFunction, 5);
                        }
                    }

                    if (endBin > startBin) {
                        scores.add(new NamedScore(s, e, v, probeName));
                    } else {
                        if (!accumulator.hasData()) accumulatedStart = s;
                        accumulatedEnd = e;
                        accumulator.add(e - s, v, probeName);
                    }

                    lastEndBin = endBin;
                }

                // Cleanup
                if (accumulator.hasData()) {
                    scores.add(getCompositeScore(accumulator, accumulatedStart, accumulatedEnd));
                }

                tile.addAllScores(scores);
            }

        }


        return tile;
    }

    private LocusScore getCompositeScore(Accumulator accumulator, int accumulatedStart, int accumulatedEnd) {
        LocusScore ls;
        if (accumulator.getNpts() == 1) {
            ls = new NamedScore(accumulatedStart, accumulatedEnd, accumulator.getRepData()[0], accumulator.getRepProbes()[0]);
        } else {
            float value = accumulator.getValue();
            ls = new CompositeScore(accumulatedStart, accumulatedEnd, value, accumulator.getRepData(),
                    accumulator.getRepProbes(), windowFunction);
        }
        return ls;

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

}
