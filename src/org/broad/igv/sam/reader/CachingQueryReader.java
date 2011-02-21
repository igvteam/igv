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
package org.broad.igv.sam.reader;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.SequenceManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentBlock;
import org.broad.igv.sam.EmptyAlignmentIterator;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LRUCache;
import org.broad.igv.util.RuntimeUtils;

import javax.swing.*;
import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.*;

/**
 * A wrapper for an AlignmentQueryReader that caches query results
 *
 * @author jrobinso
 */
public class CachingQueryReader {

    private static Logger log = Logger.getLogger(CachingQueryReader.class);

    private static Set<WeakReference<CachingQueryReader>> activeReaders = Collections.synchronizedSet(new HashSet());

      private static void cancelReaders() {
          for (WeakReference<CachingQueryReader> readerRef : activeReaders) {
              CachingQueryReader reader = readerRef.get();
              if (reader != null) {
                  reader.cancel = true;
              }
          }
          log.debug("Readers canceled");
          activeReaders.clear();
      }



    //private static final int LOW_MEMORY_THRESHOLD = 150000000;
    private static final int KB = 1000;
    private static final int MITOCHONDRIA_TILE_SIZE = 1000;
    private static int DEFAULT_TILE_SIZE = 16 * KB;
    private static int MAX_TILE_COUNT = 4;

    private String cachedChr = "";
    private int tileSize = DEFAULT_TILE_SIZE;
    private AlignmentQueryReader reader;

    private LRUCache<Integer, AlignmentTile> cache;

    private boolean cancel = false;

    public CachingQueryReader(AlignmentQueryReader reader) {
        this.reader = reader;
        cache = new LRUCache(this, MAX_TILE_COUNT);
        float fvw = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        tileSize = Math.min(DEFAULT_TILE_SIZE, (int) (fvw * KB));
    }

    public AlignmentQueryReader getWrappedReader() {
        return reader;
    }

    public void close() throws IOException {
        reader.close();
    }

    public Set<String> getSequenceNames() {
        return reader.getSequenceNames();
    }

    public SAMFileHeader getHeader() throws IOException {
        return reader.getHeader();
    }

    public CloseableIterator<Alignment> iterator() {
        return reader.iterator();
    }

    public boolean hasIndex() {
        return reader.hasIndex();
    }

    public CloseableIterator<Alignment> query(String sequence, int start, int end, List<AlignmentCounts> counts) {

        // Get the tiles covering this interval
        int startTile = (start + 1) / getTileSize(sequence);
        int endTile = end / getTileSize(sequence);    // <= inclusive
        List<AlignmentTile> tiles = getTiles(sequence, startTile, endTile);
        if (tiles.size() == 0) {
            return EmptyAlignmentIterator.getInstance();
        }

        // Count total # of records
        int recordCount = tiles.get(0).getOverlappingRecords().size();
        for (AlignmentTile t : tiles) {
            recordCount += t.getContainedRecords().size();
        }

        List<Alignment> alignments = new ArrayList(recordCount);
        alignments.addAll(tiles.get(0).getOverlappingRecords());
        for (AlignmentTile t : tiles) {
            alignments.addAll(t.getContainedRecords());
            counts.add(t.getCounts());
        }
        return new TiledIterator(start, end, alignments);
    }

    public List<AlignmentTile> getTiles(String seq, int startTile, int endTile) {

        if (!seq.equals(cachedChr)) {
            cache.clear();
            cachedChr = seq;
        }

        List<AlignmentTile> tiles = new ArrayList(endTile - startTile + 1);
        List<AlignmentTile> tilesToLoad = new ArrayList(endTile - startTile + 1);

        int tileSize = getTileSize(seq);
        for (int t = startTile; t <= endTile; t++) {
            AlignmentTile tile = cache.get(t);

            if (tile == null) {
                int start = t * tileSize;
                int end = start + tileSize;
                tile = new AlignmentTile(seq, t, start, end);
            }

            tiles.add(tile);

            // The current tile is loaded,  load any preceding tiles we have pending and clear "to load" list
            if (tile.isLoaded()) {
                if (tilesToLoad.size() > 0) {
                    boolean success = loadTiles(seq, tilesToLoad);
                    if(!success) {
                        // Loading was canceled, return what we have
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

    /**
     * Load alignments for the list of tiles
     *
     * @param chr
     * @param tiles
     * @return true if successful,  false if canceled.
     */
    private boolean loadTiles(String chr, List<AlignmentTile> tiles) {

        assert (tiles.size() > 0);

        boolean filterFailedReads = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_FILTER_FAILED_READS);
        ReadGroupFilter filter = ReadGroupFilter.getFilter();
        boolean showDuplicates = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_DUPLICATES);
        int qualityThreshold = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_QUALITY_THRESHOLD);

        //maxReadCount = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_MAX_READS);

        if (log.isDebugEnabled()) {
            int first = tiles.get(0).getTileNumber();
            int end = tiles.get(tiles.size() - 1).getTileNumber();
            log.debug("Loading tiles: " + first + "-" + end);
        }

        int start = tiles.get(0).start;
        int end = tiles.get(tiles.size() - 1).end;


        CloseableIterator<Alignment> iter = null;

        //log.debug("Loading : " + start + " - " + end);
        int alignmentCount = 0;
        try {
            //activeReaders.add(this);
            iter = reader.query(chr, start, end, false);

            int tileSize = getTileSize(chr);
            while (iter != null && iter.hasNext()) {

                if (cancel) {
                    return false;
                }

                Alignment record = iter.next();


                if (!record.isMapped() ||
                        (!showDuplicates && record.isDuplicate()) ||
                        (filterFailedReads && record.isVendorFailedRead()) ||
                        record.getMappingQuality() < qualityThreshold ||
                        (filter != null && filter.filterAlignment(record))) {
                    continue;
                }

                // Range of tile indeces that this alignment contributes to.
                int aStart = record.getAlignmentStart();
                int aEnd = record.getEnd();
                int idx0 = Math.max(0, (aStart - start) / tileSize);
                int idx1 = Math.min(tiles.size() - 1, (aEnd - start) / tileSize);

                // Loop over tiles this read overlaps
                for (int i = idx0; i <= idx1; i++) {
                    AlignmentTile t = tiles.get(i);
                    t.addRecord(record);
                }

                alignmentCount++;
                if (alignmentCount % 1000 == 0) {
                    if (cancel) return false;
                    IGVMainFrame.getInstance().setStatusBarMessage("Reads loaded: " + alignmentCount);
                    if (checkMemory() == false) {
                        cancelReaders();
                        return false;        // <=  TODO need to cancel all readers
                    }
                }
            }

            for (AlignmentTile t : tiles) {
                t.setLoaded(true);
                cache.put(t.getTileNumber(), t);
            }

            return true;

        } catch (Exception e) {
            log.error("Error loading alignment data", e);
            throw new DataLoadException("", "Error: " + e.toString());
        }

        finally {
            // reset cancel flag.  It doesn't matter how we got here,  the read is complete and this flag is reset
            // for the next time
            cancel = false;
            // activeReaders.remove(this);
            if (iter != null) {
                iter.close();
            }
            IGVMainFrame.getInstance().resetStatusMessage();
        }
    }

    private static synchronized int confirmMaxReadCount(String msg) {
        log.debug("Enter max read count");
        return JOptionPane.showConfirmDialog(null, msg, msg, JOptionPane.YES_NO_OPTION);
    }

    private static synchronized boolean checkMemory() {
        if (RuntimeUtils.getAvailableMemoryFraction() < 0.2) {
            LRUCache.clearCaches();
            System.gc();
            if (RuntimeUtils.getAvailableMemoryFraction() < 0.2) {
                String msg = "Memory is low, reading terminating.";
                MessageUtils.showMessage(msg);
                return false;
            }

        }
        return true;
    }

    /**
     * @return the tileSize
     */
    public int getTileSize(String chr) {
        if (chr.equals("M") || chr.equals("chrM") || chr.equals("MT") || chr.equals("chrMT")) {
            return MITOCHONDRIA_TILE_SIZE;
        } else {
            return tileSize;
        }
    }

    public void clearCache() {
        if (cache != null) cache.clear();
    }


    public class TiledIterator implements CloseableIterator<Alignment> {

        Iterator<Alignment> currentSamIterator;
        int end;
        Alignment nextRecord;
        int start;
        List<Alignment> alignments;

        TiledIterator(int start, int end, List<Alignment> alignments) {
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

        public Alignment next() {
            Alignment ret = nextRecord;

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
                if (nextRecord.getAlignmentStart() > end) {
                    nextRecord = null;
                }
            } else {
                nextRecord = null;
            }
        }


    }

    /**
     * Caches alignments and counts for the coverage plot.
     * <p/>
     * Notes:
     * A "bucket" is a virtual container holding all the alignments with identical start position.  The concept
     * is introduced to control the # of alignments we hold in memory for deep coverage regions.  In practice,
     * little or no information is added by displaying more than ~50X coverage.  For an average alignment length L
     * and coverage deptgh D we do not need to store more than D/L alignments at any given start position.
     */

    public static class AlignmentTile {

        private boolean loaded = false;
        private List<Alignment> containedRecords;
        private int end;
        private List<Alignment> overlappingRecords;
        private int start;
        private int tileNumber;
        private AlignmentCounts counts;

        /**
         * Maximum # of alignments to load with a given start position
         */
        int maxBucketSize;
        private int lastStart = -1;
        private int bucketCount = 0;
        private List<Alignment> currentBucket;


        AlignmentTile(String chr, int tileNumber, int start, int end) {
            this.tileNumber = tileNumber;
            this.start = start;
            this.end = end;
            containedRecords = new ArrayList(16000);
            overlappingRecords = new ArrayList();
            this.counts = new AlignmentCounts(chr, start, end);

            int maxDepth = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_MAX_LEVELS);
            // TODO -- compute this value from the data
            maxBucketSize = (maxDepth / 10) + 1;
            currentBucket = new ArrayList(5 * maxBucketSize);
        }

        public int getTileNumber() {
            return tileNumber;
        }


        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }


        public void addRecord(Alignment record) {
            if (lastStart != record.getStart()) {
                emptyBucket();
                lastStart = record.getStart();

            }
            currentBucket.add(record);
            counts.incCounts(record);
        }

        private void emptyBucket() {

            List<Alignment> sampledRecords = sampleCurrentBucket();
            for (Alignment alignment : sampledRecords) {
                int aStart = alignment.getAlignmentStart();
                int aEnd = alignment.getEnd();
                if ((aStart >= start) && (aStart < end)) {
                    containedRecords.add(alignment);
                } else if ((aEnd >= start) && (aStart < start)) {
                    overlappingRecords.add(alignment);
                }
            }
            currentBucket.clear();
        }

        private List<Alignment> sampleCurrentBucket() {

            if (currentBucket.size() < maxBucketSize) {
                return currentBucket;
            }

            Random rand = new Random(System.currentTimeMillis()); // would make this static to the class

            List subsetList = new ArrayList(maxBucketSize);
            for (int i = 0; i < maxBucketSize; i++) {
                // be sure to use Vector.remove() or you may get the same item twice
                subsetList.add(currentBucket.remove(rand.nextInt(currentBucket.size())));
            }

            return subsetList;

        }


        public List<Alignment> getContainedRecords() {
            return containedRecords;
        }


        public List<Alignment> getOverlappingRecords() {
            return overlappingRecords;
        }

        public boolean isLoaded() {
            return loaded;
        }

        public void setLoaded(boolean loaded) {
            this.loaded = loaded;

            // Empty any remaining alignments in the current bucket
            emptyBucket();

        }

        public AlignmentCounts getCounts() {
            return counts;
        }
    }

    public static class AlignmentCounts {

        String genomeId;
        //String chr;
        int start;
        int end;
        byte[] reference;
        // counts
        int[] posA;
        int[] posT;
        int[] posC;
        int[] posG;
        int[] posN;
        int[] negA;
        int[] negT;
        int[] negC;
        int[] negG;
        int[] negN;
        int[] qA;
        int[] qT;
        int[] qC;
        int[] qG;
        int[] qN;
        int[] posTotal;
        int[] negTotal;
        private int[] totalQ;
        private int maxCount = 0;

        public AlignmentCounts(String chr, int start, int end) {

            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            this.genomeId = genome.getId();
            String chrAlias = genome.getChromosomeAlias(chr);
            this.start = start;
            this.end = end;
            reference = SequenceManager.readSequence(this.genomeId, chrAlias, start, end);

            int nPts = end - start;
            posA = new int[nPts];
            posT = new int[nPts];
            posC = new int[nPts];
            posG = new int[nPts];
            posN = new int[nPts];
            posTotal = new int[nPts];
            negA = new int[nPts];
            negT = new int[nPts];
            negC = new int[nPts];
            negG = new int[nPts];
            negN = new int[nPts];
            negTotal = new int[nPts];
            qA = new int[nPts];
            qT = new int[nPts];
            qC = new int[nPts];
            qG = new int[nPts];
            qN = new int[nPts];
            totalQ = new int[nPts];
        }

        public int getTotalCount(int pos) {
            int offset = pos - start;
            if (offset < 0 || offset >= posA.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                return posTotal[offset] + negTotal[offset];

            }
        }

        public int getNegTotal(int pos) {
            int offset = pos - start;
            if (offset < 0 || offset >= posA.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                return negTotal[offset];

            }
        }

        public int getPosTotal(int pos) {
            int offset = pos - start;
            if (offset < 0 || offset >= posA.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                return posTotal[offset];

            }
        }

        public int getTotalQuality(int pos) {
            int offset = pos - start;
            if (offset < 0 || offset >= posA.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                return totalQ[offset];

            }
        }

        public int getAvgQuality(int pos) {
            int count = getTotalCount(pos);
            return count == 0 ? 0 : getTotalQuality(pos) / count;
        }

        public byte getReference(int pos) {
            if (reference == null) {
                return 0;
            }
            int offset = pos - start;
            if (offset < 0 || offset >= reference.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                return reference[offset];
            }
        }

        public int getCount(int pos, byte b) {
            int offset = pos - start;
            if (offset < 0 || offset >= posA.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                switch (b) {
                    case 'a':
                    case 'A':
                        return posA[offset] + negA[offset];
                    case 't':
                    case 'T':
                        return posT[offset] + negT[offset];
                    case 'c':
                    case 'C':
                        return posC[offset] + negC[offset];
                    case 'g':
                    case 'G':
                        return posG[offset] + negG[offset];
                    case 'n':
                    case 'N':
                        return posN[offset] + negN[offset];
                }
                log.error("Unknown nucleotide: " + b);
                return 0;
            }
        }

        public int getNegCount(int pos, byte b) {
            int offset = pos - start;
            if (offset < 0 || offset >= posA.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                switch (b) {
                    case 'a':
                    case 'A':
                        return negA[offset];
                    case 't':
                    case 'T':
                        return negT[offset];
                    case 'c':
                    case 'C':
                        return negC[offset];
                    case 'g':
                    case 'G':
                        return negG[offset];
                    case 'n':
                    case 'N':
                        return negN[offset];
                }
                log.error("Unknown nucleotide: " + b);
                return 0;
            }
        }

        public int getPosCount(int pos, byte b) {
            int offset = pos - start;
            if (offset < 0 || offset >= posA.length) {
                if (log.isDebugEnabled()) {
                    log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                }
                return 0;
            } else {
                switch (b) {
                    case 'a':
                    case 'A':
                        return posA[offset];
                    case 't':
                    case 'T':
                        return posT[offset];
                    case 'c':
                    case 'C':
                        return posC[offset];
                    case 'g':
                    case 'G':
                        return posG[offset];
                    case 'n':
                    case 'N':
                        return posN[offset];
                }
                log.error("Unknown nucleotide: " + b);
                return 0;
            }
        }

        public int getQuality(int pos, byte b) {
            int offset = pos - start;
            if (offset < 0 || offset >= posA.length) {
                log.error("Position out of range: " + pos + " (valid range - " + start + "-" + end);
                return 0;
            } else {
                switch (b) {
                    case 'a':
                    case 'A':
                        return qA[offset];
                    case 't':
                    case 'T':
                        return qT[offset];
                    case 'c':
                    case 'C':
                        return qC[offset];
                    case 'g':
                    case 'G':
                        return qG[offset];
                    case 'n':
                    case 'N':
                        return qN[offset];
                }
                log.error("Unknown nucleotide: " + posN);
                return 0;
            }
        }

        public int getAvgQuality(int pos, byte b) {
            int count = getCount(pos, b);
            return count == 0 ? 0 : getQuality(pos, b) / count;
        }


        // For alignments without blocks -- TODO refactor, this is ugly

        void incCounts(Alignment alignment) {
            int start = alignment.getAlignmentStart();
            int end = alignment.getAlignmentEnd();

            AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
            if (blocks != null) {
                for (AlignmentBlock b : blocks) {
                    // Don't count softclips
                    if (!b.isSoftClipped())
                        incCounts(b, alignment.isNegativeStrand());
                }
            } else {
                for (int pos = start; pos < end; pos++) {
                    byte q = 0;
                    incCount(pos, (byte) 'n', q, alignment.isNegativeStrand());
                }
            }
        }

        private void incCounts(AlignmentBlock block, boolean isNegativeStrand) {
            int start = block.getStart();
            byte[] bases = block.getBases();
            if (bases != null) {
                for (int i = 0; i < bases.length; i++) {
                    int pos = start + i;
                    byte q = block.getQuality(i);
                    // TODO -- handle "="
                    byte n = bases[i];
                    incCount(pos, n, q, isNegativeStrand);
                }
            }
        }

        private void incCount(int pos, byte b, byte q, boolean isNegativeStrand) {

            int offset = pos - start;
            if (offset > 0 && offset < posA.length) {
                switch (b) {
                    case 'a':
                    case 'A':
                        if (isNegativeStrand) {
                            negA[offset] = negA[offset] + 1;
                        } else {
                            posA[offset] = posA[offset] + 1;
                        }
                        qA[offset] = qA[offset] + q;
                        break;
                    case 't':
                    case 'T':
                        if (isNegativeStrand) {
                            negT[offset] = negT[offset] + 1;
                        } else {
                            posT[offset] = posT[offset] + 1;
                        }
                        qT[offset] = qT[offset] + q;
                        break;
                    case 'c':
                    case 'C':
                        if (isNegativeStrand) {
                            negC[offset] = negC[offset] + 1;
                        } else {
                            posC[offset] = posC[offset] + 1;
                        }
                        qC[offset] = qC[offset] + q;
                        break;
                    case 'g':
                    case 'G':
                        if (isNegativeStrand) {
                            negG[offset] = negG[offset] + 1;
                        } else {
                            posG[offset] = posG[offset] + 1;
                        }
                        qG[offset] = qG[offset] + q;
                        break;
                    case 'n':
                    case 'N':
                        if (isNegativeStrand) {
                            negN[offset] = negN[offset] + 1;
                        } else {
                            posN[offset] = posN[offset] + 1;
                        }
                        qN[offset] = qN[offset] + q;

                }
                if (isNegativeStrand) {
                    negTotal[offset] = negTotal[offset] + 1;
                } else {
                    posTotal[offset] = posTotal[offset] + 1;
                }
                totalQ[offset] = totalQ[offset] + q;

                maxCount = Math.max(posTotal[offset] + negTotal[offset], maxCount);
            }
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
         * @return the totalQ
         */
        public int[] getTotalQ() {
            return totalQ;
        }

        /**
         * @return the maxCount
         */
        public int getMaxCount() {
            return maxCount;
        }
    }
}


