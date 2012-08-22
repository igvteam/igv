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


import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.ReadGroupFilter;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ObjectCache;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.collections.LRUCache;

import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.*;

/**
 * A wrapper for an AlignmentQueryReader that caches query results
 *
 * @author jrobinso
 */
public class NoncachingQueryReader {

    private static Logger log = Logger.getLogger(CachingQueryReader.class);

    //private static final int LOW_MEMORY_THRESHOLD = 150000000;
    private static final int KB = 1000;
    private static final int MITOCHONDRIA_TILE_SIZE = 1000;
    private static int MAX_TILE_COUNT = 10;
    private static Set<WeakReference<NoncachingQueryReader>> activeReaders = Collections.synchronizedSet(new HashSet());

    /**
     * Flag to mark a corrupt index.  Without this attempted reads will continue in an infinite loop
     */
    private boolean corruptIndex = false;

    private float visibilityWindow = 16;    // Visibility window,  in KB
    private String cachedChr = "";
    private int tileSize;
    private AlignmentReader reader;
    private boolean cancel = false;
    private boolean pairedEnd = false;


    // Map of read group -> paired end stats

    //private PairedEndStats peStats;

    private static void cancelReaders() {
        for (WeakReference<NoncachingQueryReader> readerRef : activeReaders) {
            NoncachingQueryReader reader = readerRef.get();
            if (reader != null) {
                reader.cancel = true;
            }
        }
        log.debug("Readers canceled");
        activeReaders.clear();
    }


    public NoncachingQueryReader(AlignmentReader reader) {
        this.reader = reader;
        activeReaders.add(new WeakReference<NoncachingQueryReader>(this));
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

    /*
               List<AlignmentCounts> counts = new ArrayList();
                List<DownsampledInterval> downsampledIntervals = new ArrayList<DownsampledInterval>();
                List<SpliceJunctionFeature> spliceJunctions = null;
                if (showSpliceJunctions) {
                    spliceJunctions = new ArrayList<SpliceJunctionFeature>();
                }

                Iterator<Alignment> iter = reader.query(sequence, intervalStart, intervalEnd, counts,
                        spliceJunctions, downsampledIntervals, downsampleOptions, peStats, bisulfiteContext);

                 final AlignmentPacker alignmentPacker = new AlignmentPacker();

                LinkedHashMap<String, List<AlignmentInterval.Row>> alignmentRows = alignmentPacker.packAlignments(iter,
                        intervalEnd, renderOptions);

                AlignmentInterval loadedInterval = new AlignmentInterval(chr, intervalStart, intervalEnd,
                        alignmentRows, counts, spliceJunctions, downsampledIntervals, renderOptions);

     */

    public AlignmentInterval loadInterval(String chr, int start, int end, boolean showSpliceJunctions,
                                          AlignmentTrack.RenderOptions renderOptions,
                                          Map<String, PEStats> peStats,
                                          AlignmentTrack.BisulfiteContext bisulfiteContext) {


        AlignmentDataManager.DownsampleOptions downsampleOptions = new AlignmentDataManager.DownsampleOptions();

        AlignmentTile t = loadTiles(chr, start, end, downsampleOptions, peStats, bisulfiteContext);

        List<Alignment> alignments =  t.getAlignments();

        List<SpliceJunctionFeature> spliceJunctions = t.getSpliceJunctionFeatures();

        List<AlignmentCounts> counts = new ArrayList();
        counts.add(t.getCounts());

        List<DownsampledInterval> downsampledIntervals = t.getDownsampledIntervals();

        // Since we (potentially) downsampled,  we need to sort
        Collections.sort(alignments, new AlignmentSorter());

        Iterator<Alignment> iter =  alignments.iterator();

        final AlignmentPacker alignmentPacker = new AlignmentPacker();

        LinkedHashMap<String, List<AlignmentInterval.Row>> alignmentRows = alignmentPacker.packAlignments(iter,
                end, renderOptions);

        return new AlignmentInterval(chr, start, end, alignmentRows, counts, spliceJunctions, downsampledIntervals,
                renderOptions);

    }

    private Iterator<Alignment> query(String chr, int start, int end,
                                     List<AlignmentCounts> counts,
                                     List<SpliceJunctionFeature> spliceJunctionFeatures,
                                     List<DownsampledInterval> downsampledIntervals,
                                     AlignmentDataManager.DownsampleOptions downsampleOptions,
                                     Map<String, PEStats> peStats,
                                     AlignmentTrack.BisulfiteContext bisulfiteContext) {


        AlignmentTile t = loadTiles(chr, start, end, downsampleOptions, peStats, bisulfiteContext);

        if (t == null) {
            return (new ArrayList<Alignment>()).iterator();
        }

        List<Alignment> alignments =  t.getAlignments();

        if (spliceJunctionFeatures != null) {
            List<SpliceJunctionFeature> tmp = t.getSpliceJunctionFeatures();
            if (tmp != null) spliceJunctionFeatures.addAll(tmp);
        }


        counts.add(t.getCounts());
        downsampledIntervals.addAll(t.getDownsampledIntervals());


        // Since we added in 2 passes, and downsampled,  we need to sort
        Collections.sort(alignments, new AlignmentSorter());

        return alignments.iterator();
    }


    private AlignmentTile loadTiles(String chr, int start, int end,
                                    AlignmentDataManager.DownsampleOptions downsampleOptions,
                                    Map<String, PEStats> peStats,
                                    AlignmentTrack.BisulfiteContext bisulfiteContext) {

        AlignmentTile t = new AlignmentTile(start, end, downsampleOptions, bisulfiteContext);


        //assert (tiles.size() > 0);
        if (corruptIndex) {
            return t;
        }

        boolean filterFailedReads = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_FILTER_FAILED_READS);
        ReadGroupFilter filter = ReadGroupFilter.getFilter();
        boolean showDuplicates = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_DUPLICATES);
        int qualityThreshold = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_QUALITY_THRESHOLD);

        CloseableIterator<Alignment> iter = null;

        //log.debug("Loading : " + start + " - " + end);
        int alignmentCount = 0;
        WeakReference<NoncachingQueryReader> ref = new WeakReference(this);
        try {
            ObjectCache<String, Alignment> mappedMates = new ObjectCache<String, Alignment>(1000);
            ObjectCache<String, Alignment> unmappedMates = new ObjectCache<String, Alignment>(1000);


            activeReaders.add(ref);
            iter = reader.query(chr, start, end, false);

            while (iter != null && iter.hasNext()) {

                if (cancel) {
                    return t;
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

                t.addRecord(record);

                alignmentCount++;
                int interval = Globals.isTesting() ? 100000 : 1000;
                if (alignmentCount % interval == 0) {
                    if (cancel) return null;
                    MessageUtils.setStatusBarMessage("Reads loaded: " + alignmentCount);
                    if (checkMemory() == false) {
                        cancelReaders();
                        return t;        // <=  TODO need to cancel all readers
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
            t.setLoaded(true);

            return t;

        } catch (java.nio.BufferUnderflowException e) {
            // This almost always indicates a corrupt BAM index, or less frequently a corrupt bam file
            corruptIndex = true;
            MessageUtils.showMessage("<html>Error encountered querying alignments: " + e.toString() +
                    "<br>This is often caused by a corrupt index file.");
            return null;

        } catch (Exception e) {
            log.error("Error loading alignment data", e);
            MessageUtils.showMessage("<html>Error encountered querying alignments: " + e.toString());
            return null;
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
     * Does this file contain paired end data?  Assume not until proven otherwise.
     */
    public boolean isPairedEnd() {
        return pairedEnd;
    }

    public Set<String> getPlatforms() {
        return reader.getPlatforms();
    }

    /**
     * Caches alignments, coverage, splice junctions, and downsampled intervals
     */

    public static class AlignmentTile {

        private boolean loaded = false;
        private int end;
        private int start;
        private AlignmentCounts counts;
        private List<Alignment> alignments;
        private List<DownsampledInterval> downsampledIntervals;
        private List<SpliceJunctionFeature> spliceJunctionFeatures;
        private SpliceJunctionHelper spliceJunctionHelper;

        private boolean downsample;
        private int samplingWindowSize;
        private int samplingDepth;
        private SamplingBucket currentSamplingBucket;

        private static final Random RAND = new Random(System.currentTimeMillis());


        AlignmentTile(int start, int end,
                      AlignmentDataManager.DownsampleOptions downsampleOptions,
                      AlignmentTrack.BisulfiteContext bisulfiteContext) {
            this.start = start;
            this.end = end;
            alignments = new ArrayList(16000);
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
                spliceJunctionFeatures = new ArrayList<SpliceJunctionFeature>(100);
                spliceJunctionHelper = new SpliceJunctionHelper();
            }

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
                alignments.add(alignment);
        }

        public List<Alignment> getAlignments() {
            return alignments;
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
                        spliceJunctionFeatures.add(f);
                }
            }
            spliceJunctionHelper = null;
        }


        public List<SpliceJunctionFeature> getSpliceJunctionFeatures() {
            return spliceJunctionFeatures;
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


}


