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

package org.broad.igv.methyl;

import org.apache.log4j.Logger;
import org.broad.igv.util.collections.LRUCache;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author jrobinso
 * @date Jun 24, 2010
 */
public class CachingMethylSource implements MethylDataSource {

    private static Logger log = Logger.getLogger(CachingMethylSource.class);
    private static int DEFAULT_TILE_COUNT = 4;
    private int binSize;

    MethylDataSource reader;
    LRUCache<String, Bin> cache;


    public CachingMethylSource(MethylDataSource reader, int binSize) {
        this(reader, DEFAULT_TILE_COUNT, binSize);
    }


    public CachingMethylSource(MethylDataSource reader, int tileCount, int binSize) {
        this.reader = reader;
        this.cache = new LRUCache(tileCount);
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


    public Iterator<MethylScore> query(String chr, int start, int end) {

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
        int recordCount = 0;
        for (Bin t : tiles) {
            recordCount += t.getContainedRecords().size();
        }

        List<MethylScore> alignments = new ArrayList(recordCount);
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
        Iterator<MethylScore> iter = null;

        //log.debug("Loading : " + start + " - " + end);
        int featureCount = 0;
        long t0 = System.currentTimeMillis();
        try {


            iter = reader.query(seq, start, end);

            while (iter != null && iter.hasNext()) {
                MethylScore record = iter.next();

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

                    // A bin size == 0 means use a single bin for the entire chromosome.  This is a confusing convention.
                    if (binSize == 0 || ((aStart >= t.start) && (aStart < t.end))) {
                        t.containedRecords.add(record);
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

        } finally {
            if (iter != null) {
                //iter.close();
            }
            //IGV.getInstance().resetStatusMessage();
        }
    }


    static class Bin {

        private boolean loaded = false;
        private int start;
        private int end;
        private int binNumber;
        private List<MethylScore> containedRecords;

        Bin(int binNumber, int start, int end) {
            this.binNumber = binNumber;
            this.start = start;
            this.end = end;
            containedRecords = new ArrayList(1000);
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

        public List<MethylScore> getContainedRecords() {
            return containedRecords;
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
    public class BinIterator implements Iterator<MethylScore> {

        Iterator<MethylScore> currentSamIterator;
        int end;
        MethylScore nextRecord;
        int start;
        List<MethylScore> alignments;

        BinIterator(int start, int end, List<MethylScore> alignments) {
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

        public MethylScore next() {
            MethylScore ret = nextRecord;

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

