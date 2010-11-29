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


package org.broad.igv.track.tribble;

import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureSource;
import org.broad.tribble.iterators.CloseableTribbleIterator;
import org.broad.tribble.util.LRUCache;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * @author jrobinso
 * @date Jun 24, 2010
 */
public class CachingFeatureReader implements org.broad.tribble.FeatureSource {

    private static Logger log = Logger.getLogger(CachingFeatureReader.class);
    private static int maxBinCount = 10;
    private static int defaultBinSize = 16000; // <= 16 kb

    private int binSize;
    FeatureSource reader;
    LRUCache<String, Bin> cache;


    public CachingFeatureReader(FeatureSource reader) {
        this(reader, maxBinCount, defaultBinSize);
    }


    public CachingFeatureReader(FeatureSource reader, int binCount, int binSize) {
        this.reader = reader;
        this.cache = new LRUCache(binCount);
        this.binSize = binSize;
    }


    /**
     * Set the bin size.   This invalidates the cache.
     *
     * @param newSize
     */
    public void setBinSize(int newSize) {
        this.binSize = newSize;
        cache.clear();

    }


    public List<String> getSequenceNames() {
        return reader.getSequenceNames();
    }

    public Object getHeader() {
        return reader.getHeader();
    }

    /**
     * Return an iterator over the entire file.  Nothing to cache,  just delegate to the wrapped reader
     *
     * @throws java.io.IOException
     */
    public CloseableTribbleIterator iterator() throws IOException {
        return reader.iterator();
    }

    public void close() throws IOException {
        cache.clear();
        reader.close();
    }


    public CloseableTribbleIterator query(String chr, int start, int end) throws IOException {

        // A binSize of zero => use a single bin for the entire chromosome
        int startBin = 0;
        int endBin = 0;    // <= inclusive
        if (binSize > 0) {
            startBin = start / binSize;
            endBin = end / binSize;    // <= inclusive
        }
        List<Bin> tiles = getBins(chr, startBin, endBin);

        if (tiles.size() == 0) {
            return null;
        }

        // Count total # of records
        int recordCount = tiles.get(0).getOverlappingRecords().size();
        for (Bin t : tiles) {
            recordCount += t.getContainedRecords().size();
        }

        List<Feature> alignments = new ArrayList(recordCount);
        alignments.addAll(tiles.get(0).getOverlappingRecords());
        for (Bin t : tiles) {
            alignments.addAll(t.getContainedRecords());
        }
        return new BinIterator(start, end, alignments);
    }


    /**
     * Return loaded tiles that span the query interval.
     *
     * @param seq
     * @param startBin
     * @param endBin
     * @return
     */
    private List<Bin> getBins(String seq, int startBin, int endBin) {

        List<Bin> tiles = new ArrayList(endBin - startBin + 1);
        List<Bin> tilesToLoad = new ArrayList(endBin - startBin + 1);

        for (int t = startBin; t <= endBin; t++) {
            String key = seq + "_" + t;
            Bin tile = cache.get(key);

            if (tile == null) {
                if (log.isDebugEnabled()) {
                    log.debug("Tile cache miss: " + t);
                }
                int start = t * binSize;
                int end = start + binSize;
                tile = new Bin(t, start, end);
                cache.put(key, tile);
            }

            tiles.add(tile);

            // The current tile is loaded,  load any preceding tiles we have pending
            if (tile.isLoaded()) {
                if (tilesToLoad.size() > 0) {
                    if (!loadTiles(seq, tilesToLoad)) {
                        return tiles;
                    }
                }
                tilesToLoad.clear();
            } else {
                tilesToLoad.add(tile);
            }
        }

        if (tilesToLoad.size() > 0) {
            loadTiles(seq, tilesToLoad);
        }

        return tiles;
    }

    private boolean loadTiles(String seq, List<Bin> tiles) {

        assert (tiles.size() > 0);

        if (log.isDebugEnabled()) {
            int first = tiles.get(0).getBinNumber();
            int end = tiles.get(tiles.size() - 1).getBinNumber();
            log.debug("Loading tiles: " + first + "-" + end);
        }

        // Convert start to 1-based coordinates
        int start = tiles.get(0).start + 1;
        int end = tiles.get(tiles.size() - 1).end;
        Iterator<Feature> iter = null;

        //log.debug("Loading : " + start + " - " + end);
        int featureCount = 0;
        long t0 = System.currentTimeMillis();
        try {
            iter = reader.query(seq, start, end);

            while (iter != null && iter.hasNext()) {
                Feature record = iter.next();

                // Range of tile indeces that this alignment contributes to.
                int aStart = record.getStart();
                int aEnd = record.getEnd();
                int idx0 = 0;
                int idx1 = 0;
                if (binSize > 0) {
                    idx0 = Math.max(0, (aStart - start) / binSize);
                    idx1 = Math.min(tiles.size() - 1, (aEnd - start) / binSize);
                }

                // Loop over tiles this read overlaps
                for (int i = idx0; i <= idx1; i++) {
                    Bin t = tiles.get(i);

                    if ((aStart >= t.start) && (aStart < t.end)) {
                        t.containedRecords.add(record);
                    } else if ((aEnd >= t.start) && (aStart < t.start)) {
                        t.overlappingRecords.add(record);
                    }
                }
            }

            for (Bin t : tiles) {
                t.setLoaded(true);
            }
            if (log.isDebugEnabled()) {
                long dt = System.currentTimeMillis() - t0;
                long rate = dt == 0 ? Long.MAX_VALUE : featureCount / dt;
                log.debug("Loaded " + featureCount + " reads in " + dt + "ms.  (" + rate + " reads/ms)");
            }
            return true;

        } catch (IOException e) {
            log.error("IOError loading alignment data", e);

            // TODO -- do something about this,  how do we want to handle this exception?
            throw new RuntimeException(e);
        }

        finally {
            if (iter != null) {
                //iter.close();
            }
            //IGVMainFrame.getInstance().resetStatusMessage();
        }
    }


    static class Bin {

        private boolean loaded = false;
        private int start;
        private int end;
        private int binNumber;
        private List<Feature> containedRecords;
        private List<Feature> overlappingRecords;

        Bin(int binNumber, int start, int end) {
            this.binNumber = binNumber;
            this.start = start;
            this.end = end;
            containedRecords = new ArrayList(1000);
            overlappingRecords = new ArrayList(100);
        }

        public int getBinNumber() {
            return binNumber;
        }


        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }

        public List<Feature> getContainedRecords() {
            return containedRecords;
        }

        public List<Feature> getOverlappingRecords() {
            return overlappingRecords;
        }

        public boolean isLoaded() {
            return loaded;
        }

        public void setLoaded(boolean loaded) {
            this.loaded = loaded;
        }

    }

    /**
     * TODO -- this is a pointeless class.  It would make sense if it actually took tiles, instead of the collection
     * TODO -- of alignments.
     */
    public class BinIterator implements CloseableTribbleIterator {

        Iterator<Feature> currentSamIterator;
        int end;
        Feature nextRecord;
        int start;
        List<Feature> alignments;

        BinIterator(int start, int end, List<Feature> alignments) {
            this.alignments = alignments;
            this.start = start;
            this.end = end;
            currentSamIterator = alignments.iterator();
            advanceToFirstRecord();
        }

        public void close() {
            // No-op
        }

        public boolean hasNext() {
            return nextRecord != null;
        }

        public Feature next() {
            Feature ret = nextRecord;

            advanceToNextRecord();

            return ret;
        }

        public void remove() {
            // ignored
        }

        private void advanceToFirstRecord() {
            advanceToNextRecord();
        }

        private void advanceToNextRecord() {
            advance();

            while ((nextRecord != null) && (nextRecord.getEnd() < start)) {
                advance();
            }
        }

        private void advance() {
            if (currentSamIterator.hasNext()) {
                nextRecord = currentSamIterator.next();
                if (nextRecord.getStart() > end) {
                    nextRecord = null;
                }
            } else {
                nextRecord = null;
            }
        }

        public Iterator iterator() {
            return this;
        }
    }
}

