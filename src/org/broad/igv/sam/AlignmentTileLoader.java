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

package org.broad.igv.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.ReadGroupFilter;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.ObjectCache;
import org.broad.igv.util.RuntimeUtils;

import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.*;

/**
 * A wrapper for an AlignmentQueryReader that caches query results
 *
 * @author jrobinso
 */
public class AlignmentTileLoader {

    private static Logger log = Logger.getLogger(AlignmentTileLoader.class);

    private static Set<WeakReference<AlignmentTileLoader>> activeLoaders = Collections.synchronizedSet(new HashSet());

    /**
     * Flag to mark a corrupt index.  Without this attempted reads will continue in an infinite loop
     */
    private boolean corruptIndex = false;

    private AlignmentReader reader;
    private boolean cancel = false;
    private boolean pairedEnd = false;

    static void cancelReaders() {
        for (WeakReference<AlignmentTileLoader> readerRef : activeLoaders) {
            AlignmentTileLoader reader = readerRef.get();
            if (reader != null) {
                reader.cancel = true;
            }
        }
        log.debug("Readers canceled");
        activeLoaders.clear();
    }


    public AlignmentTileLoader(AlignmentReader reader) {
        this.reader = reader;
        activeLoaders.add(new WeakReference<AlignmentTileLoader>(this));
    }

    public void close() throws IOException {
        reader.close();
    }

    public SAMFileHeader getFileHeader() {
        return this.reader.getFileHeader();
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


    AlignmentTile loadTile(String chr, int start, int end,
                           SpliceJunctionHelper spliceJunctionHelper,
                           AlignmentDataManager.DownsampleOptions downsampleOptions,
                           Map<String, PEStats> peStats,
                           AlignmentTrack.BisulfiteContext bisulfiteContext,
                           ProgressMonitor monitor) {

        AlignmentTile t = new AlignmentTile(start, end, spliceJunctionHelper, downsampleOptions, bisulfiteContext);


        //assert (tiles.size() > 0);
        if (corruptIndex) {
            return t;
        }

        final PreferenceManager prefMgr = PreferenceManager.getInstance();
        boolean filterFailedReads = prefMgr.getAsBoolean(PreferenceManager.SAM_FILTER_FAILED_READS);
        boolean filterSecondaryAlignments = prefMgr.getAsBoolean(PreferenceManager.SAM_FILTER_SECONDARY_ALIGNMENTS);
        boolean filterSupplementaryAlignments = prefMgr.getAsBoolean(PreferenceManager.SAM_FILTER_SUPPLEMENTARY_ALIGNMENTS);
        ReadGroupFilter filter = ReadGroupFilter.getFilter();
        boolean showDuplicates = prefMgr.getAsBoolean(PreferenceManager.SAM_SHOW_DUPLICATES);
        int qualityThreshold = prefMgr.getAsInt(PreferenceManager.SAM_QUALITY_THRESHOLD);

        CloseableIterator<Alignment> iter = null;

        //log.debug("Loading : " + start + " - " + end);
        int alignmentCount = 0;
        WeakReference<AlignmentTileLoader> ref = new WeakReference(this);
        try {
            ObjectCache<String, Alignment> mappedMates = new ObjectCache<String, Alignment>(1000);
            ObjectCache<String, Alignment> unmappedMates = new ObjectCache<String, Alignment>(1000);


            activeLoaders.add(ref);
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
                        (filterSecondaryAlignments && !record.isPrimary()) ||
                        (filterSupplementaryAlignments && record.isSupplementary()) ||
                        record.getMappingQuality() < qualityThreshold ||
                        (filter != null && filter.filterAlignment(record))) {
                    continue;
                }

                t.addRecord(record);

                alignmentCount++;
                int interval = Globals.isTesting() ? 100000 : 1000;
                if (alignmentCount % interval == 0) {
                    if (cancel) return null;
                    String msg = "Reads loaded: " + alignmentCount;
                    MessageUtils.setStatusBarMessage(msg);
                    if (monitor != null) {
                        monitor.updateStatus(msg);
                    }
                    if (memoryTooLow()) {
                        if (monitor != null) monitor.fireProgressChange(100);
                        cancelReaders();
                        t.finish();
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
                double minPercentile = prefMgr.getAsFloat(PreferenceManager.SAM_MIN_INSERT_SIZE_PERCENTILE);
                double maxPercentile = prefMgr.getAsFloat(PreferenceManager.SAM_MAX_INSERT_SIZE_PERCENTILE);
                for (PEStats stats : peStats.values()) {
                    stats.compute(minPercentile, maxPercentile);
                }
            }

            // Clean up any remaining unmapped mate sequences
            for (String mappedMateName : mappedMates.getKeys()) {
                Alignment mappedMate = mappedMates.get(mappedMateName);
                if (mappedMate != null) {
                    Alignment mate = unmappedMates.get(mappedMate.getReadName());
                    if (mate != null) {
                        mappedMate.setMateSequence(mate.getReadSequence());
                    }
                }
            }
            t.finish();


        } catch (java.nio.BufferUnderflowException e) {
            // This almost always indicates a corrupt BAM index, or less frequently a corrupt bam file
            corruptIndex = true;
            MessageUtils.showMessage("<html>Error encountered querying alignments: " + e.toString() +
                    "<br>This is often caused by a corrupt index file.");

        } catch (Exception e) {
            log.error("Error loading alignment data", e);
            MessageUtils.showMessage("<html>Error encountered querying alignments: " + e.toString());
        } finally {
            // reset cancel flag.  It doesn't matter how we got here,  the read is complete and this flag is reset
            // for the next time
            cancel = false;
            activeLoaders.remove(ref);

            if (monitor != null) {
                monitor.fireProgressChange(100);
            }

            if (iter != null) {
                iter.close();
            }
            if (!Globals.isHeadless()) {
                IGV.getInstance().resetStatusMessage();
            }
        }

        return t;

    }


    private static synchronized boolean memoryTooLow() {
        if (RuntimeUtils.getAvailableMemoryFraction() < 0.2) {
            System.gc();
            if (RuntimeUtils.getAvailableMemoryFraction() < 0.2) {
                String msg = "Memory is low, reading terminating.";
                MessageUtils.showMessage(msg);
                return true;
            }

        }
        return false;
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
        private SpliceJunctionHelper spliceJunctionHelper;

        private static final Random RAND = new Random();

        private boolean downsample;
        private int samplingWindowSize;
        private int samplingDepth;

        private int currentSamplingWindowStart = -1;
        private int curEffSamplingWindowDepth = 0;
        //Although not strictly necessary, we keep track of the currentDownsampledInterval
        //for easy incrementing
        private DownsampledInterval currentDownsampledInterval;

        /**
         * We keep a data structure of alignments which can efficiently
         * look up by string (read name) or index (for random replacement)
         */
        IndexableMap<String, Alignment> imAlignments;

        private int downsampledCount = 0;
        private int offset = 0;

        AlignmentTile(int start, int end,
                      SpliceJunctionHelper spliceJunctionHelper,
                      AlignmentDataManager.DownsampleOptions downsampleOptions,
                      AlignmentTrack.BisulfiteContext bisulfiteContext) {
            this.start = start;
            this.end = end;
            this.downsampledIntervals = new ArrayList<DownsampledInterval>();
            long seed = System.currentTimeMillis();
            //System.out.println("seed: " + seed);
            RAND.setSeed(seed);

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

            this.spliceJunctionHelper = spliceJunctionHelper;

            if (this.downsample) {
                imAlignments = new IndexableMap<String, Alignment>(8000);
            } else {
                alignments = new ArrayList<Alignment>(16000);
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
        // * @param alignment
         */
        public void addRecord(Alignment alignment) {

           // alignment = new ReducedMemoryAlignment(alignment);

            counts.incCounts(alignment);

            if (spliceJunctionHelper != null) {
                spliceJunctionHelper.addAlignment(alignment);
            }

            if (downsample) {
                final int alignmentStart = alignment.getAlignmentStart();
                int currentSamplingBucketEnd = currentSamplingWindowStart + samplingWindowSize;
                if (currentSamplingWindowStart < 0 || alignmentStart >= currentSamplingBucketEnd) {
                    setCurrentSamplingBucket(alignmentStart);
                }

                attemptAddRecordDownsampled(alignment);

            } else {
                alignments.add(alignment);
            }

            alignment.finish();
        }

        /**
         * Attempt to add this alignment. The alignment is definitely added if there is another
         * read with the same name. Typically this other read is a mate pair, but it could also be a secondary alignment
         * or chimeric read.
         * If we haven't seen another read like that, the record is added with some probability according to
         * reservoir sampling
         *
         * @param alignment
         */
        private void attemptAddRecordDownsampled(Alignment alignment) {
            String readName = alignment.getReadName();
            //A simple way to turn off the same-readName-checking is to replace the read name with a random string
            //so that there are no repeats
            //readName = String.format("%s%d", readName, RAND.nextInt());

            //There are 3 possibilities: other-kept, other-rejected, other-unknown (haven't seen)
            //If we kept or rejected the another read with the same name, we do the same for this one
            boolean hasRead = imAlignments.containsKey(readName);
            if (hasRead) {
                List<Alignment> mateAlignments = imAlignments.get(readName);
                boolean haveOther = mateAlignments != null;
                if (haveOther) {
                    //We keep the alignment if others have been kept
                    imAlignments.append(readName, alignment);
                } else {
                    currentDownsampledInterval.incCount();
                }
            } else {
                if (curEffSamplingWindowDepth < samplingDepth) {
                    imAlignments.append(readName, alignment);
                    curEffSamplingWindowDepth++;
                } else {
                    double samplingProb = ((double) samplingDepth) / (samplingDepth + downsampledCount + 1);
                    if (RAND.nextDouble() < samplingProb) {
                        int rndInt = (int) (RAND.nextDouble() * (samplingDepth - 1));
                        int idx = offset + rndInt;
                        // Replace random record with this one
                        List<Alignment> removedValues = imAlignments.replace(idx, readName, alignment);
                        incrementDownsampledIntervals(removedValues);
                    } else {
                        //Mark that record was not kept
                        imAlignments.markNull(readName);
                        currentDownsampledInterval.incCount();
                    }
                    downsampledCount++;
                }
            }
        }

        private void setCurrentSamplingBucket(int alignmentStart) {
            curEffSamplingWindowDepth = 0;
            downsampledCount = 0;
            currentSamplingWindowStart = alignmentStart;
            offset = imAlignments.size();

            int currentSamplingBucketEnd = currentSamplingWindowStart + samplingWindowSize;
            currentDownsampledInterval = new DownsampledInterval(alignmentStart, currentSamplingBucketEnd, 0);
            downsampledIntervals.add(currentDownsampledInterval);
        }

        private void incrementDownsampledIntervals(List<Alignment> removedValues) {
            if (removedValues == null) return;
            for (Alignment al : removedValues) {
                DownsampledInterval interval = findDownsampledInterval(al);
                if (interval != null) interval.incCount();
            }
        }

        private DownsampledInterval findDownsampledInterval(Alignment al) {
            return findDownsampledInterval(al, 0, downsampledIntervals.size());
        }

        /**
         * Attempt to find appropriate DownsampledInterval by recursive binary search
         *
         * @param al
         * @param startInd Start search index, 0-based inclusive
         * @param endInd   End search index, 0-based exclusive
         * @return
         */
        private DownsampledInterval findDownsampledInterval(Alignment al, int startInd, int endInd) {
            //Length-0 search space
            if (startInd == endInd) {
                return null;
            }
            final int midInd = (startInd + endInd) / 2;
            final DownsampledInterval curInterval = downsampledIntervals.get(midInd);
            if (al.getStart() >= curInterval.getStart() && al.getStart() < curInterval.getEnd()) {
                //Found
                return curInterval;
            }

            if (al.getStart() >= curInterval.getEnd()) {
                startInd = midInd + 1;
            } else {
                endInd = midInd;
            }
            return findDownsampledInterval(al, startInd, endInd);
        }

        /**
         * Sort the alignments by start position, and filter {@code downsampledIntervals}.
         * This will have the same results as if no downsampling occurred, although will incur
         * extra computational cost
         */
        private void sortFilterDownsampled() {
            if ((this.alignments == null || this.alignments.size() == 0) && this.downsample) {
                this.alignments = imAlignments.getAllValues();
                imAlignments.clear();
            }

            Comparator<Alignment> alignmentSorter = new Comparator<Alignment>() {
                public int compare(Alignment alignment, Alignment alignment1) {
                    return alignment.getStart() - alignment1.getStart();
                }
            };
            Collections.sort(this.alignments, alignmentSorter);

            //Only keep the intervals for which count > 0
            List<DownsampledInterval> tmp = new ArrayList<DownsampledInterval>(this.downsampledIntervals.size());
            for (DownsampledInterval interval : this.downsampledIntervals) {
                if (interval.getCount() > 0) {
                    tmp.add(interval);
                }
            }
            this.downsampledIntervals = tmp;
        }

        public List<Alignment> getAlignments() {

            if(alignments == null) {
                finish();   // TODO -- I'm not sure this should ever happen
            }
            return alignments;
        }

        public List<DownsampledInterval> getDownsampledIntervals() {
            return downsampledIntervals;
        }

        public void finish() {
            //If we downsampled,  we need to sort
            if (downsample) {
                sortFilterDownsampled();
            }
            finalizeSpliceJunctions();
            counts.finish();
        }

        public AlignmentCounts getCounts() {
            return counts;
        }


        private void finalizeSpliceJunctions() {
            if (spliceJunctionHelper != null) {
           //     spliceJunctionHelper.finish();
            }
        }


        /**
         * Map-like structure designed to be accessible both by key, and by numeric index
         * Multiple values are stored for each key, and a list is returned
         * If the value for a key is set as null, nothing can be added
         * <p/>
         * Intended to support downsampling, where if a read name is added and then removed
         * we don't want to add the read pair
         *
         * @param <K>
         * @param <V>
         */
        private class IndexableMap<K, V> {
            private HashMap<K, List<V>> map;
            private List<K> list;

            IndexableMap(int size) {
                this.map = new HashMap<K, List<V>>(size);
                this.list = new ArrayList<K>(size);
            }

            public List<V> get(K key) {
                return map.get(key);
            }

            /**
             * Append a value for the specified key, unless
             * the current value is null. If the current value is
             * null, it's a no-op.
             *
             * @param key
             * @param value
             * @return Whether the element was added
             */
            public boolean append(K key, V value) {
                if (!map.containsKey(key)) {
                    addNewValueToMap(key, value);
                    return list.add(key);
                } else {
                    List<V> curList = map.get(key);
                    if (curList == null) return false;
                    return curList.add(value);
                }
            }

            public List<V> markNull(K key) {
                return map.put(key, null);
            }

            private void addNewValueToMap(K key, V value) {
                List<V> curList = new ArrayList<V>(2);
                curList.add(value);
                map.put(key, curList);
            }

            /**
             * Place the specified {@code key} and {@code value} in the map,
             * at index {@code index}.
             * <p/>
             * In the unlikely event that {@code key} is already
             * at {@code index}, {@code value} will be appended
             *
             * @param index
             * @param key
             * @param value
             * @return Whether the replacement actually happened
             */
            public List<V> replace(int index, K key, V value) {
                checkSize(index);
                K oldKey = list.get(index);
                if (!oldKey.equals(key)) {
                    //Remove the old key from map, and make sure nothing else gets put there
                    List<V> oldValue = markNull(oldKey);
                    addNewValueToMap(key, value);
                    list.set(index, key);
                    return oldValue;
                } else {
                    append(key, value);
                    return null;
                }
            }

            public int size() {
                return list.size();
            }

            private void checkSize(int index) {
                if (index >= size()) {
                    throw new IllegalArgumentException("index " + index + " greater than current size" + size());
                }
            }

            public List<V> getAllValues() {
                List<V> allValues = new ArrayList<V>(2 * size());
                for (K k : list) {
                    allValues.addAll(map.get(k));
                }
                return allValues;
            }

            public boolean containsKey(K key) {
                return map.containsKey(key);
            }

            public void clear() {
                map.clear();
                list.clear();
            }
        }

    }

}


