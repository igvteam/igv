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

package org.broad.igv.sam;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.ReadGroupFilter;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LRUCache;
import org.broad.igv.util.ObjectCache;
import org.broad.igv.util.RuntimeUtils;
import org.broad.tribble.Feature;

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

    //private static final int LOW_MEMORY_THRESHOLD = 150000000;
    private static final int KB = 1000;
    private static final int MITOCHONDRIA_TILE_SIZE = 1000;
    private static int MAX_TILE_COUNT = 10;
    private static Set<WeakReference<CachingQueryReader>> activeReaders = Collections.synchronizedSet(new HashSet());

    /**
     * Flag to mark a corrupt index.  Without this attempted reads will continue in an infinite loop
     */
    private boolean corruptIndex = false;

    private float visibilityWindow = 16;    // Visibility window,  in KB
    private String cachedChr = "";
    private int tileSize;
    private AlignmentReader reader;
    private boolean cancel = false;
    private LRUCache<Integer, AlignmentTile> cache;
    private boolean pairedEnd = false;


    // Map of read group -> paired end stats

    //private PairedEndStats peStats;

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


    public CachingQueryReader(AlignmentReader reader) {
        this.reader = reader;
        activeReaders.add(new WeakReference<CachingQueryReader>(this));
        updateCache();
    }

    public static void visibilityWindowChanged() {
        for (WeakReference<CachingQueryReader> readerRef : activeReaders) {
            CachingQueryReader reader = readerRef.get();
            if (reader != null) {
                reader.updateCache();
            }
        }
    }

    private void updateCache() {
        float fvw = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);

        // If the visibility window has changed by more than a factor of 2 change cache tile size
        float ratio = fvw / visibilityWindow;
        if (cache == null || (ratio < 0.5 || ratio > 2)) {
            // Set tile size to  the visibility window
            tileSize = (int) (fvw * KB);
            cache = new LRUCache(this, MAX_TILE_COUNT);
            visibilityWindow = fvw;
        }
    }

    public AlignmentReader getWrappedReader() {
        return reader;
    }

    public void close() throws IOException {
        reader.close();
    }

    public List<String> getSequenceNames() {
        return reader.getSequenceNames();
    }

    public CloseableIterator<Alignment> iterator() {
        return reader.iterator();
    }

    public boolean hasIndex() {
        return reader.hasIndex();
    }

    public CloseableIterator<Alignment> query(String sequence, int start, int end,
                                              List<AlignmentCounts> counts,
                                              List<SpliceJunctionFeature> spliceJunctionFeatures,
                                              List<DownsampledInterval> downsampledIntervals,
                                              AlignmentDataManager.DownsampleOptions downsampleOptions,
                                              Map<String, PEStats> peStats,
                                              AlignmentTrack.BisulfiteContext bisulfiteContext) {

        // Get the tiles covering this interval
        int startTile = (start + 1) / getTileSize(sequence);
        int endTile = end / getTileSize(sequence);    // <= inclusive

        List<AlignmentTile> tiles = getTiles(sequence, startTile, endTile, downsampleOptions, peStats, bisulfiteContext);
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

        if (spliceJunctionFeatures != null) {
            List<SpliceJunctionFeature> tmp = tiles.get(0).getOverlappingSpliceJunctionFeatures();
            if (tmp != null) spliceJunctionFeatures.addAll(tmp);
        }

        for (AlignmentTile t : tiles) {
            alignments.addAll(t.getContainedRecords());
            counts.add(t.getCounts());
            downsampledIntervals.addAll(t.getDownsampledIntervals());

            if (spliceJunctionFeatures != null) {
                List<SpliceJunctionFeature> tmp = t.getContainedSpliceJunctionFeatures();
                if (tmp != null) spliceJunctionFeatures.addAll(tmp);
            }
        }

        // Since we added in 2 passes, and downsampled,  we need to sort
        Collections.sort(alignments, new AlignmentSorter());

        return new TiledIterator(start, end, alignments);
    }

    public List<AlignmentTile> getTiles(String seq, int startTile, int endTile,
                                        AlignmentDataManager.DownsampleOptions downsampleOptions,
                                        Map<String, PEStats> peStats,
                                        AlignmentTrack.BisulfiteContext bisulfiteContext) {

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

                tile = new AlignmentTile(t, start, end, downsampleOptions, bisulfiteContext);
            }

            tiles.add(tile);

            // The current tile is loaded,  load any preceding tiles we have pending and clear "to load" list
            if (tile.isLoaded()) {
                if (tilesToLoad.size() > 0) {
                    boolean success = loadTiles(seq, tilesToLoad, peStats);
                    if (!success) {
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
            loadTiles(seq, tilesToLoad, peStats);
        }

        return tiles;
    }

    /**
     * Load alignments for the list of tiles
     *
     * @param chr     Only tiles on this chromosome will be loaded
     * @param tiles
     * @param peStats
     * @return true if successful,  false if canceled.
     */
    private boolean loadTiles(String chr, List<AlignmentTile> tiles, Map<String, PEStats> peStats) {

        //assert (tiles.size() > 0);
        if (corruptIndex) {
            return false;
        }

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
        WeakReference<CachingQueryReader> ref = new WeakReference(this);
        try {
            ObjectCache<String, Alignment> mappedMates = new ObjectCache<String, Alignment>(1000);
            ObjectCache<String, Alignment> unmappedMates = new ObjectCache<String, Alignment>(1000);


            activeReaders.add(ref);
            iter = reader.query(chr, start, end, false);

            int tileSize = getTileSize(chr);
            while (iter != null && iter.hasNext()) {

                if (cancel) {
                    return false;
                }

                Alignment record = iter.next();

                // Set mate sequence of unmapped mates
                // Put a limit on the total size of this collection.
                String readName = record.getReadName();
                if (record.isPaired()) {
                    pairedEnd = true;
                    if (record.isMapped()) {
                        if (!record.getMate().isMapped()) {
                            // record is mapped, mate is not
                            Alignment mate = unmappedMates.get(readName);
                            if (mate == null) {
                                mappedMates.put(readName, record);
                            } else {
                                record.setMateSequence(mate.getReadSequence());
                                unmappedMates.remove(readName);
                                mappedMates.remove(readName);
                            }

                        }
                    } else if (record.getMate().isMapped()) {
                        // record not mapped, mate is
                        Alignment mappedMate = mappedMates.get(readName);
                        if (mappedMate == null) {
                            unmappedMates.put(readName, record);
                        } else {
                            mappedMate.setMateSequence(record.getReadSequence());
                            unmappedMates.remove(readName);
                            mappedMates.remove(readName);
                        }
                    }
                }


                if (!record.isMapped() || (!showDuplicates && record.isDuplicate()) ||
                        (filterFailedReads && record.isVendorFailedRead()) ||
                        record.getMappingQuality() < qualityThreshold ||
                        (filter != null && filter.filterAlignment(record))) {
                    continue;
                }

                // Range of tile indices that this alignment contributes to.
                int aStart = record.getStart();
                int aEnd = record.getEnd();
                int idx0 = Math.max(0, (aStart - start) / tileSize);
                int idx1 = Math.min(tiles.size() - 1, (aEnd - start) / tileSize);

                // Loop over tiles this read overlaps
                for (int i = idx0; i <= idx1; i++) {
                    AlignmentTile t = null;
                    t = tiles.get(i);
                    t.addRecord(record);
                }

                alignmentCount++;
                int interval = Globals.isTesting() ? 100000 : 1000;
                if (alignmentCount % interval == 0) {
                    if (cancel) return false;
                    MessageUtils.setStatusBarMessage("Reads loaded: " + alignmentCount);
                    if (checkMemory() == false) {
                        cancelReaders();
                        return false;        // <=  TODO need to cancel all readers
                    }
                }

                // Update pe stats
                if (peStats != null && record.isPaired() && record.isProperPair()) {
                    String lb = record.getLibrary();
                    if (lb == null) lb = "null";
                    PEStats stats = peStats.get(lb);
                    if (stats == null) {
                        stats = new PEStats(lb);
                        peStats.put(lb, stats);
                    }
                    stats.update(record);

                }
            }
            // End iteration over alignments

            // Compute peStats
            if (peStats != null) {
                // TODO -- something smarter re the percentiles.  For small samples these will revert to min and max
                double minPercentile = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MIN_INSERT_SIZE_PERCENTILE);
                double maxPercentile = PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_MAX_INSERT_SIZE_PERCENTILE);
                for (PEStats stats : peStats.values()) {
                    stats.compute(minPercentile, maxPercentile);
                }
            }

            // Clean up any remaining unmapped mate sequences
            for (String mappedMateName : mappedMates.getKeys()) {
                Alignment mappedMate = mappedMates.get(mappedMateName);
                Alignment mate = unmappedMates.get(mappedMate.getReadName());
                if (mate != null) {
                    mappedMate.setMateSequence(mate.getReadSequence());
                }
            }
            mappedMates = null;
            unmappedMates = null;

            for (AlignmentTile t : tiles) {
                t.setLoaded(true);
                cache.put(t.getTileNumber(), t);
            }

            return true;

        } catch (java.nio.BufferUnderflowException e) {
            // This almost always indicates a corrupt BAM index, or less frequently a corrupt bam file
            corruptIndex = true;
            MessageUtils.showMessage("<html>Error encountered querying alignments: " + e.toString() +
                    "<br>This is often caused by a corrupt index file.");
            return false;

        } catch (Exception e) {
            log.error("Error loading alignment data", e);
            MessageUtils.showMessage("<html>Error encountered querying alignments: " + e.toString());
            return false;
        } finally {
            // reset cancel flag.  It doesn't matter how we got here,  the read is complete and this flag is reset
            // for the next time
            cancel = false;
            activeReaders.remove(ref);
            if (iter != null) {
                iter.close();
            }
            if (!Globals.isHeadless()) {
                IGV.getInstance().resetStatusMessage();
            }
        }
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
     * @param chr Chromosome name
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

    /**
     * Does this file contain paired end data?  Assume not until proven otherwise.
     */
    public boolean isPairedEnd() {
        return pairedEnd;
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

            //We use exclusive end
            while ((nextRecord != null) && (nextRecord.getEnd() <= start)) {
                advance();
            }
        }

        private void advance() {
            if (currentSamIterator.hasNext()) {
                nextRecord = currentSamIterator.next();

                if (nextRecord.getStart() >= end) {
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
     * and coverage depth D we do not need to store more than D/L alignments at any given start position.
     */

    //

    public static class AlignmentTile {

        private boolean loaded = false;
        private int end;
        private int start;
        private int tileNumber;
        private AlignmentCounts counts;
        private List<Alignment> containedRecords;
        private List<Alignment> overlappingRecords;
        private List<DownsampledInterval> downsampledIntervals;
        private List<SpliceJunctionFeature> containedSpliceJunctionFeatures;
        private List<SpliceJunctionFeature> overlappingSpliceJunctionFeatures;
        private SpliceJunctionHelper spliceJunctionHelper;

        private boolean downsample;
        private int samplingWindowSize;
        private int samplingDepth;
        private SamplingBucket currentSamplingBucket;

        private static final Random RAND = new Random(System.currentTimeMillis());


        AlignmentTile(int tileNumber, int start, int end,
                      AlignmentDataManager.DownsampleOptions downsampleOptions,
                      AlignmentTrack.BisulfiteContext bisulfiteContext) {
            this.tileNumber = tileNumber;
            this.start = start;
            this.end = end;
            containedRecords = new ArrayList(16000);
            overlappingRecords = new ArrayList();
            downsampledIntervals = new ArrayList();

            // Use a sparse array for large regions  (> 10 mb)
            if ((end - start) > 10000000) {
                this.counts = new SparseAlignmentCounts(start, end, bisulfiteContext);
            } else {
                this.counts = new DenseAlignmentCounts(start, end, bisulfiteContext);
            }

            // Set the max depth, and the max depth of the sampling bucket.
            if (downsampleOptions == null) {
                // Use default settings (from preferences)
                downsampleOptions = new AlignmentDataManager.DownsampleOptions();
            }
            this.downsample = downsampleOptions.isDownsample();
            this.samplingWindowSize = downsampleOptions.getSampleWindowSize();
            this.samplingDepth = Math.max(1, downsampleOptions.getMaxReadCount());

            // TODO -- only if splice junctions are on
            if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_JUNCTION_TRACK)) {
                containedSpliceJunctionFeatures = new ArrayList<SpliceJunctionFeature>(100);
                overlappingSpliceJunctionFeatures = new ArrayList<SpliceJunctionFeature>(100);
                spliceJunctionHelper = new SpliceJunctionHelper();
            }

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

        int ignoredCount = 0;    // <= just for debugging

        /**
         * Add an alignment record to this tile.  This record is not necessarily retained after down-sampling.
         *
         * @param alignment
         */
        public void addRecord(Alignment alignment) {

            counts.incCounts(alignment);

            if (downsample) {
                final int alignmentStart = alignment.getAlignmentStart();
                if (currentSamplingBucket == null || alignmentStart >= currentSamplingBucket.end) {
                    if (currentSamplingBucket != null) {
                        emptyBucket();
                    }
                    int end = alignmentStart + samplingWindowSize;
                    currentSamplingBucket = new SamplingBucket(alignmentStart, end);
                }

                if (spliceJunctionHelper != null) {
                    spliceJunctionHelper.addAlignment(alignment);
                }

                currentSamplingBucket.add(alignment);
            } else {
                allocateAlignment(alignment);
            }
        }

        private void emptyBucket() {
            if (currentSamplingBucket == null) {
                return;
            }
            //List<Alignment> sampledRecords = sampleCurrentBucket();
            for (Alignment alignment : currentSamplingBucket.getAlignments()) {
                allocateAlignment(alignment);
            }

            if (currentSamplingBucket.isSampled()) {
                DownsampledInterval interval = new DownsampledInterval(currentSamplingBucket.start,
                        currentSamplingBucket.end, currentSamplingBucket.downsampledCount);
                downsampledIntervals.add(interval);
            }

            currentSamplingBucket = null;

        }

        private void allocateAlignment(Alignment alignment) {
            int aStart = alignment.getStart();
            int aEnd = alignment.getEnd();
            if ((aStart >= start) && (aStart < end)) {
                containedRecords.add(alignment);
            } else if ((aEnd > start) && (aStart < start)) {
                overlappingRecords.add(alignment);
            }
        }

        public List<Alignment> getContainedRecords() {
            return containedRecords;
        }

        public List<Alignment> getOverlappingRecords() {
            return overlappingRecords;
        }

        public List<DownsampledInterval> getDownsampledIntervals() {
            return downsampledIntervals;
        }

        public boolean isLoaded() {
            return loaded;
        }

        public void setLoaded(boolean loaded) {
            this.loaded = loaded;

            if (loaded) {
                // Empty any remaining alignments in the current bucket
                emptyBucket();
                currentSamplingBucket = null;
                finalizeSpliceJunctions();
                counts.finish();
            }
        }

        public AlignmentCounts getCounts() {
            return counts;
        }


        private void finalizeSpliceJunctions() {
            if (spliceJunctionHelper != null) {
                spliceJunctionHelper.finish();
                List<SpliceJunctionFeature> features = spliceJunctionHelper.getFeatures();
                for (SpliceJunctionFeature f : features) {
                    if (f.getStart() >= start) {
                        containedSpliceJunctionFeatures.add(f);
                    } else {
                        overlappingSpliceJunctionFeatures.add(f);
                    }
                }
            }
            spliceJunctionHelper = null;
        }


        public List<SpliceJunctionFeature> getContainedSpliceJunctionFeatures() {
            return containedSpliceJunctionFeatures;
        }

        public List<SpliceJunctionFeature> getOverlappingSpliceJunctionFeatures() {
            return overlappingSpliceJunctionFeatures;
        }


        private class SamplingBucket {
            int start;
            int end;
            int downsampledCount = 0;
            List<Alignment> alignments;


            private SamplingBucket(int start, int end) {
                this.start = start;
                this.end = end;
                alignments = new ArrayList(samplingDepth);
            }


            public void add(Alignment alignment) {
                // If the current bucket is < max depth we keep it.  Otherwise,  keep with probability == samplingProb
                // If we have the mate in the bucket already, always keep it.
                if (alignments.size() < samplingDepth) {
                    alignments.add(alignment);
                } else {
                    double samplingProb = ((double) samplingDepth) / (samplingDepth + downsampledCount + 1);
                    if (RAND.nextDouble() < samplingProb) {
                        int idx = (int) (RAND.nextDouble() * (alignments.size() - 1));
                        // Replace random record with this one
                        alignments.set(idx, alignment);
                    }
                    downsampledCount++;

                }

            }

            public List<Alignment> getAlignments() {
                return alignments;
            }

            public boolean isSampled() {
                return downsampledCount > 0;
            }

            public int getDownsampledCount() {
                return downsampledCount;
            }
        }
    }

    private static class AlignmentSorter implements Comparator<Alignment> {
        public int compare(Alignment alignment, Alignment alignment1) {
            return alignment.getStart() - alignment1.getStart();
        }
    }


    public static class DownsampledInterval implements Feature {
        private int start;
        private int end;
        private int count;

        public DownsampledInterval(int start, int end, int count) {
            this.start = start;
            this.end = end;
            this.count = count;
        }

        public String toString() {
            return start + "-" + end + " (" + count + ")";
        }

        public int getCount() {
            return count;
        }

        public int getEnd() {
            return end;
        }

        public int getStart() {
            return start;
        }

        public String getChr() {
            return null;
        }

        public String getValueString() {
            return "Interval [" + start + "-" + end + "] <br>" + count + " reads removed.";
        }
    }


}


