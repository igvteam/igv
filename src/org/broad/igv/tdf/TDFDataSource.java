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
package org.broad.igv.tdf;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.LRUCache;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * @author jrobinso
 */
public class TDFDataSource implements DataSource {

    private static Logger log = Logger.getLogger(TDFDataSource.class);
    // TODO -- read from file
    int maxPrecomputedZoom = 6;
    TDFReader reader;
    private int trackNumber = 0;
    String trackName;
    LRUCache<String, List<LocusScore>> summaryScoreCache = new LRUCache(this, 20);
    Genome genome;
    Interval currentInterval;
    WindowFunction windowFunction = WindowFunction.mean;
    List<WindowFunction> availableFunctions;

    private boolean aggregateLikeBins = true;

    float normalizationFactor = 1.0f;


    public TDFDataSource(TDFReader reader, int trackNumber, String trackName) {

        boolean normalizeCounts = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.NORMALIZE_COVERAGE);

        // TODO -- a single reader will be shared across data sources
        this.trackNumber = trackNumber;
        this.trackName = trackName;
        this.reader = reader;
        this.availableFunctions = reader.getWindowFunctions();

        TDFGroup rootGroup = reader.getGroup("/");
        try {
            maxPrecomputedZoom = Integer.parseInt(rootGroup.getAttribute("maxZoom"));
        } catch (Exception e) {
            log.error("Error reading attribute 'maxZoom'", e);
        }
        try {
            String dataGenome = rootGroup.getAttribute("genome");
            // TODO -- throw exception if data genome != current genome
            genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
        } catch (Exception e) {
            log.error("Unknown genome " + rootGroup.getAttribute("genome"));
            throw new RuntimeException("Unknown genome " + rootGroup.getAttribute("genome"));
        }
        try {
            if (normalizeCounts) {
                String totalCountString = rootGroup.getAttribute("totalCount");
                if (totalCountString != null) {
                    int totalCount = Integer.parseInt(totalCountString);
                    normalizationFactor = 1.0e6f / totalCount;
                }
            }
        } catch (Exception e) {
            log.error("Error reading attribute 'maxZoom'", e);
        }

    }

    public String getPath() {
        return reader == null ? null : reader.getPath();
    }

    public String getTrackName() {
        return trackName;
    }

    public double getDataMax() {
        return reader.getUpperLimit() * normalizationFactor;
    }

    public double getDataMin() {
        return reader.getLowerLimit() * normalizationFactor;
    }

    public void setAggregateLikeBins(boolean aggregateLikeBins) {
        this.aggregateLikeBins = aggregateLikeBins;
    }

    class Interval {

        String chr;
        private int start;
        private int end;
        private int zoom;
        private List<LocusScore> scores;

        public Interval(String chr, int start, int end, int zoom, List<LocusScore> scores) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.zoom = zoom;
            this.scores = scores;
        }

        public boolean contains(String chr, int s, int e, int zoom) {
            return chr.equals(this.chr) && zoom == this.zoom && s >= getStart() && e <= getEnd();
        }

        /**
         * @return the start
         */
        public int getStart() {
            return start;
        }

        /**
         * @return the end
         */
        public int getEnd() {
            return end;
        }

        /**
         * @return the scores
         */
        public List<LocusScore> getScores() {
            return scores;
        }
    }

    private List<LocusScore> getCachedSummaryScores(String chr, int zoom, int tileNumber, double tileWidth) {

        String key = chr + "_" + zoom + "_" + tileNumber + "_" + windowFunction;

        List<LocusScore> scores = summaryScoreCache.get(key);
        if (scores == null) {

            int startLocation = (int) (tileNumber * tileWidth);
            int endLocation = (int) ((tileNumber + 1) * tileWidth);

            scores = getSummaryScores(chr, zoom, startLocation, endLocation);

            summaryScoreCache.put(key, scores);
        }

        return scores;

    }

    private List<LocusScore> getSummaryScores(String chr, int zoom, int startLocation, int endLocation) {
        List<LocusScore> scores;
        if (zoom <= this.maxPrecomputedZoom) {
            scores = new ArrayList(1000);
            TDFDataset ds = reader.getDataset(chr, zoom, windowFunction);
            if (ds != null) {
                List<TDFTile> tiles = ds.getTiles(startLocation, endLocation);
                if (tiles.size() > 0) {
                    for (TDFTile tile : tiles) {

                        if (tile.getSize() > 0) {
                            for (int i = 0; i < tile.getSize(); i++) {
                                float v = tile.getValue(trackNumber, i);
                                if (!Float.isNaN(v)) {
                                    v *= normalizationFactor;
                                    scores.add(new BasicScore(tile.getStartPosition(i), tile.getEndPosition(i), v));
                                }
                            }
                        }
                    }
                }
            }
        } else {

            int chrLength = getChrLength(chr);
            if (chrLength == 0) {
                return Collections.emptyList();
            }
            endLocation = Math.min(endLocation, chrLength);
            // By definition there are 2^z tiles per chromosome, and 700 bins per tile, where z is the zoom level.
            int nTiles = (int) Math.pow(2, zoom);
            float binSize = (((float) chrLength) / nTiles) / 700;

            scores = computeSummaryScores(chr, startLocation, endLocation, zoom, binSize);
        }
        return scores;
    }


    public int getChrLength(String chr) {
        final Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
        if (chr.equals(Globals.CHR_ALL)) {
            return (int) (genome.getLength() / 1000);
        } else {
            Chromosome c = genome.getChromosome(chr);
            return c == null ? 0 : c.getLength();
        }
    }

    private List<LocusScore> computeSummaryScores(String chr, int startLocation, int endLocation, int zoom, float binSize) {

        List<LocusScore> scores = new ArrayList(1000);

        String dsName = "/" + chr + "/raw";

        TDFDataset rawDataset = reader.getDataset(dsName);
        if (rawDataset != null) {
            List<TDFTile> rawTiles = rawDataset.getTiles(startLocation, endLocation);
            if (rawTiles.size() > 0) {

                if (windowFunction == WindowFunction.none) {
                    for (TDFTile rawTile : rawTiles) {
                        // Tile of raw data
                        if (rawTile != null && rawTile.getSize() > 0) {

                            for (int i = 0; i < rawTile.getSize(); i++) {
                                int s = rawTile.getStartPosition(i);
                                int e = Math.max(s, rawTile.getEndPosition(i) - 1);

                                if (e < startLocation) {
                                    continue;
                                } else if (s > endLocation) {
                                    break;
                                }
                                float v = rawTile.getValue(trackNumber, i);
                                if (!Float.isNaN(v)) {
                                    v *= normalizationFactor;
                                }
                                scores.add(new BasicScore(s, e, v));
                            }
                        }
                    }


                } else {

                    // The minimum bins size is 1 bp
                    int nBins = (int) ((endLocation - startLocation) / binSize + 1);

                    Bin[] bins = new Bin[nBins];

                    for (TDFTile rawTile : rawTiles) {
                        // Tile of raw data
                        if (rawTile != null && rawTile.getSize() > 0) {

                            for (int i = 0; i < rawTile.getSize(); i++) {
                                int s = rawTile.getStartPosition(i);
                                int e = Math.max(s, rawTile.getEndPosition(i) - 1);

                                if (e < startLocation) {
                                    continue;
                                } else if (s > endLocation) {
                                    break;
                                }


                                String probeName = rawTile.getName(i);
                                float v = rawTile.getValue(trackNumber, i);
                                if (!Float.isNaN(v)) {
                                    v *= normalizationFactor;
                                    int startBin = (int) Math.max(0, ((s - startLocation) / binSize));
                                    int endBin = (int) Math.min(bins.length - 1, ((e - startLocation) / binSize));
                                    for (int b = startBin; b <= endBin; b++) {
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
                                if (aggregateLikeBins && currentBin.isExtension(bins[b])) {
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
                }
            }
        }
        return scores;
    }


    /*
     *
    Chromosome chromosome = IGV.getInstance().getGenomeManager().getGenome().getChr(chr);
    if(chromosome == null) {
    return null;
    }
    double tileWidth = chromosome.getLength() / (Math.pow(2.0, zoom));
     */

    public synchronized List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {

        Chromosome chromosome = genome.getChromosome(chr);


        // If we are in gene list view bypass caching.
        if (FrameManager.isGeneListMode()) {

            return getSummaryScores(chr, zoom, startLocation, endLocation);

        } else {

            ArrayList scores = new ArrayList();

            // TODO -- this whole section could be computed once and stored,  it is only a function of the genome, chr, and zoom level.
            double tileWidth = 0;
            if (chr.equals(Globals.CHR_ALL)) {
                tileWidth = (genome.getLength() / 1000.0);
            } else {
                if (chromosome != null) {
                    tileWidth = chromosome.getLength() / ((int) Math.pow(2.0, zoom));
                }
            }
            if (tileWidth == 0) {
                return null;
            }

            int startTile = (int) (startLocation / tileWidth);
            int endTile = (int) (endLocation / tileWidth);
            for (int t = startTile; t <= endTile; t++) {
                List<LocusScore> cachedScores = getCachedSummaryScores(chr, zoom, t, tileWidth);
                if (cachedScores != null) {
                    for (LocusScore s : cachedScores) {
                        if (s.getEnd() >= startLocation) {
                            scores.add(s);
                        } else if (s.getStart() > endLocation) {
                            break;
                        }
                    }
                }

            }

            return scores;
        }
    }


    public TrackType getTrackType() {
        return reader.getTrackType();
    }

    public void setWindowFunction(WindowFunction wf) {
        this.windowFunction = wf;
    }

    public boolean isLogNormalized() {
        //return false;
        return getDataMin() < 0;
    }

    public void refreshData(long timestamp) {
        // ignored
    }

    public WindowFunction getWindowFunction() {
        return windowFunction;
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return availableFunctions;
    }

}

