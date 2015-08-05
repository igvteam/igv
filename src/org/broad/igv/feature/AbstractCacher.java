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


package org.broad.igv.feature;

import org.apache.log4j.Logger;
import org.broad.igv.util.collections.LRUCache;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;


/**
 * Class to handle caching data from any source of features
 * Subclasses must override queryRaw
 *
 * @author jrobinso
 * @date Jun 24, 2010
 */
public abstract class AbstractCacher {

    private static Logger log = Logger.getLogger(AbstractCacher.class);

    protected int binSize = Integer.MAX_VALUE;
    protected LRUCache<String, Bin> cache;


    public AbstractCacher(int binCount, int binSize) {
        this.cache = new LRUCache(binCount);
        setBinSize(binSize);
    }

    /**
     * Obtain data from underlying source
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    protected abstract Iterator<Feature> queryRaw(String chr, int start, int end) throws IOException;



    /**
     * Set the bin size.   This invalidates the cache.
     *
     * @param newSize
     */
    public void setBinSize(int newSize) {
        this.binSize = newSize == 0 ? Integer.MAX_VALUE : newSize;  // A binSize of zero => use a single bin for the entire chromosome
        cache.clear();

    }

    public void close() throws IOException {
        cache.clear();
    }

    /**
     * Query the cached data, refreshing from raw data as necessary
     *
     * @param chr
     * @param start
     * @param end
     * @return
     * @throws IOException
     */
    public Iterator<Feature> queryCached(String chr, int start, int end) throws IOException {

        int  startBin = start / binSize;
        int endBin = end / binSize;    // <= inclusive

        List<Bin> tiles = getBins(chr, startBin, endBin);

        if (tiles.size() == 0) {
            return Collections.<Feature>emptyList().iterator();
        }

        // Count total # of records
        int recordCount = tiles.get(0).getOverlappingRecords().size();
        for (Bin t : tiles) {
            recordCount += t.getContainedRecords().size();
        }

        List<Feature> features = new ArrayList(recordCount);
        features.addAll(tiles.get(0).getOverlappingRecords());
        for (Bin t : tiles) {
            features.addAll(t.getContainedRecords());
        }
        return new BinIterator(start, end, features);
    }


    /**
     * Return loaded tiles that span the query interval.
     * <p/>
     * We synchronize this method because different threads might be using
     * the same source. Synchronizing here ensures that data
     * is loaded as few times as possible (first caller loads it into the cache,
     * the second caller accesses it from there) as well as preventing bugs stemming
     * from multiple thread access
     *
     * @param seq
     * @param startBin
     * @param endBin
     * @return
     */
    private synchronized List<Bin> getBins(String seq, int startBin, int endBin) {

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

            iter = queryRaw(seq, start, end);

            while (iter != null && iter.hasNext()) {
                Feature record = iter.next();

                // Range of tile indices that this feature contributes to.
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

                    // A bin size == 0 means use a single bin for the entire chromosome.  This is a confusing convention.
                    if (binSize == 0 || ((aStart >= t.start) && (aStart < t.end))) {
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
            log.error("IOError loading feature data", e);

            // TODO -- do something about this,  how do we want to handle this exception?
            throw new RuntimeException(e);
        } finally {
            if (iter != null) {
                //iter.close();
            }
            //IGV.getInstance().resetStatusMessage();
        }
    }


    private static class Bin {

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
     *
     */
    private class BinIterator implements CloseableTribbleIterator {

        Iterator<Feature> currentFeatureIterator;
        int end;
        Feature nextRecord;
        int start;
        List<Feature> features;

        BinIterator(int start, int end, List<Feature> features) {
            this.features = features;
            this.start = start;
            this.end = end;
            currentFeatureIterator = features.iterator();
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
            if (currentFeatureIterator.hasNext()) {
                nextRecord = currentFeatureIterator.next();
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

