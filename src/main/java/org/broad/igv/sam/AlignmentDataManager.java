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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.event.RefreshEvent;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentTrack.SortOption;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.track.RenderContext;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.net.MalformedURLException;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * Manages data loading for a single alignment file.  Shared between alignment, coverage, and junction
 * tracks.
 */

public class AlignmentDataManager implements IGVEventObserver {

    private static Logger log = Logger.getLogger(AlignmentDataManager.class);


    private AlignmentTrack alignmentTrack;
    private CoverageTrack coverageTrack;
    private Set<Track> subscribedTracks;

    private List<AlignmentInterval> intervalCache;
    private ResourceLocator locator;
    private HashMap<String, String> chrMappings = new HashMap();
    private AlignmentTileLoader reader;
    private Map<String, PEStats> peStats;
    private SpliceJunctionHelper.LoadOptions loadOptions;
    private Object loadLock = new Object();


    public AlignmentDataManager(ResourceLocator locator, Genome genome) throws IOException {
        this.locator = locator;
        // The time-gated limit for an AWS signed URL has expired, we need to re-sign the URL with the newly acquired
        // access token, otherwise we will face an Access Denied error. CheckReader() provides a very low overhead
        // mechanism to refresh expired presigned URLs.
        reader = new AlignmentTileLoader(AlignmentReaderFactory.getReader(locator));
        peStats = new HashMap();
        initLoadOptions();
        initChrMap(genome);
        intervalCache = Collections.synchronizedList(new ArrayList<>());
        subscribedTracks = Collections.synchronizedSet(new HashSet<>());

        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(RefreshEvent.class, this);
    }

    public void receiveEvent(Object event) {
        if (event instanceof FrameManager.ChangeEvent) {
            trimCache();
        } else if (event instanceof RefreshEvent) {
            clear();
        } else {
            log.info("Unknown event type: " + event.getClass());
        }
    }

    public void subscribe(Track track) {
        subscribedTracks.add(track);
    }

    public void unsubscribe(Track track) {
        subscribedTracks.remove(track);
        if (subscribedTracks.isEmpty()) {
            dumpAlignments();
            IGVEventBus.getInstance().unsubscribe(this);
        }
    }

    void initLoadOptions() {
        this.loadOptions = new SpliceJunctionHelper.LoadOptions();
    }

    /**
     * Create an alias -> chromosome lookup map.  Enables loading BAM files that use alternative names for chromosomes
     * (e.g. 1 -> chr1,  etc).
     */
    private void initChrMap(Genome genome) throws IOException {

        if (genome != null) {

            // Build a chr size -> name lookup table.   We will assume sizes are unique.  This will be used if no alias
            // is defined for a sequence.
            Map<Long, String> inverseDict = null;
            Map<String, Long> sequenceDictionary = checkReader().getSequenceDictionary();

            if (sequenceDictionary != null) {

                Set<Long> nonUnique = new HashSet<>();
                Set<Long> seen = new HashSet<>();
                // First find sequences whose size are not unique,  we'll filter these
                for(Long size : sequenceDictionary.values()) {
                    if(seen.contains(size)) {
                        nonUnique.add(size);
                    } else{
                        seen.add(size);
                    }
                }

                inverseDict = new HashMap<>();

                for (Chromosome chromosome : genome.getChromosomes()) {

                    Long size = (long) chromosome.getLength();
                    if (!nonUnique.contains(size)) {
                        if (inverseDict.containsKey(size)) {
                            inverseDict.remove(size);
                            nonUnique.add(size);
                        } else {
                            inverseDict.put(size, chromosome.getName());
                        }
                    }
                }
            }


            List<String> seqNames = checkReader().getSequenceNames();
            if (seqNames != null) {
                for (String seq : seqNames) {

                    if (genome.isKnownChr(seq)) {
                        String chr = genome.getCanonicalChrName(seq);
                        chrMappings.put(chr, seq);
                    } else if (sequenceDictionary != null) {
                        Long size = sequenceDictionary.get(seq);
                        String chr = inverseDict.get(size);
                        if (chr != null) {
                            chrMappings.put(chr, seq);
                        }
                    }
                }
            }
        }
    }

    public boolean hasMatchingSequences() {
        return chrMappings.size() > 0;
    }

    public AlignmentTileLoader getReader() {
        return checkReader();
    }

    public ResourceLocator getLocator() {
        return locator;
    }

    public Map<String, PEStats> getPEStats() {
        return peStats;
    }

    public boolean isPairedEnd() {
        return checkReader().isPairedEnd();
    }

    public boolean hasYCTags() {
        return checkReader().hasYCTags();
    }

    public boolean hasIndex() {
        return checkReader().hasIndex();
    }

    public void setAlignmentTrack(AlignmentTrack alignmentTrack) {
        this.alignmentTrack = alignmentTrack;
    }

    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public double getMinVisibleScale() {
        return getVisibilityWindow() / 700;
    }

    public double getVisibilityWindow() {
        return getPreferences().getAsFloat(SAM_MAX_VISIBLE_RANGE) * 1000;
    }

    private IGVPreferences getPreferences() {
        String category =  NULL_CATEGORY;
        AlignmentTrack.ExperimentType experimentType = getExperimentType();
        if(experimentType == AlignmentTrack.ExperimentType.RNA) {
            category = RNA;
        } else if(experimentType == AlignmentTrack.ExperimentType.THIRD_GEN) {
            category = THIRD_GEN;
        }
        return PreferencesManager.getPreferences(category);
    }


    /**
     * The set of sequences found in the file.
     * May be null
     *
     * @return
     */
    public List<String> getSequenceNames() throws IOException {
        return checkReader().getSequenceNames();
    }

    public Map<String, Long> getSequenceDictionary() {
        return checkReader().getSequenceDictionary();
    }


    public AlignmentInterval getLoadedInterval(ReferenceFrame frame) {

        for (AlignmentInterval interval : intervalCache) {
            if (interval.contains(frame.getCurrentRange())) {
                return interval;
            }
        }

        return null;
    }

    /**
     * Sort rows group by group
     *
     * @param option
     * @param location
     */
    public boolean sortRows(SortOption option, ReferenceFrame frame, double location, String tag) {

        AlignmentInterval interval = getLoadedInterval(frame);
        if (interval == null) {
            return false;
        } else {
            PackedAlignments packedAlignments = interval.getPackedAlignments();
            if (packedAlignments == null) {
                return false;
            }

            for (List<Row> alignmentRows : packedAlignments.values()) {
                for (Row row : alignmentRows) {
                    row.updateScore(option, location, interval, tag);
                }
                Collections.sort(alignmentRows);
            }
            return true;
        }
    }

    public void setViewAsPairs(boolean option, AlignmentTrack.RenderOptions renderOptions) {
        if (option == renderOptions.isViewPairs()) {
            return;
        }
        renderOptions.setViewPairs(option);
        packAlignments(renderOptions);
    }

    /**
     * Repack currently loaded alignments across frames
     * All relevant intervals must be loaded
     *
     * @param renderOptions
     * @return Whether repacking was performed
     */
    void packAlignments(AlignmentTrack.RenderOptions renderOptions) {
        for (AlignmentInterval interval : intervalCache) {
            interval.packAlignments(renderOptions);
        }
    }


    public boolean isLoaded(ReferenceFrame frame) {
        return getLoadedInterval(frame) != null;
    }

    public  void  load(ReferenceFrame frame,
                     AlignmentTrack.RenderOptions renderOptions,
                     boolean expandEnds) {

        if (frame.getChrName().equals(Globals.CHR_ALL) || frame.getScale() > getMinVisibleScale())
            return; // should not happen

        if (isLoaded(frame)) {
            return;  // Already loaded
        }
        
        synchronized (loadLock) {

            if (isLoaded(frame)) {
                return;  // Already loaded
            }

            Range range = frame.getCurrentRange();
            final String chr = frame.getChrName();
            final int start = range.getStart();
            final int end = range.getEnd();
            int adjustedStart = start;
            int adjustedEnd = end;

            // Expand the interval by the lesser of  +/- a 2 screens, or max visible range
            int windowSize = Math.min(4 * (end - start), PreferencesManager.getPreferences().getAsInt(SAM_MAX_VISIBLE_RANGE) * 1000);
            int center = (end + start) / 2;
            int expand = Math.max(end - start, windowSize / 2);

            if (expandEnds) {
                adjustedStart = Math.max(0, Math.min(start, center - expand));
                adjustedEnd = Math.max(end, center + expand);
            }

            AlignmentInterval loadedInterval = loadInterval(chr, adjustedStart, adjustedEnd, renderOptions);

            trimCache();

            intervalCache.add(loadedInterval);

            packAlignments(renderOptions);
        }

        //  IGVEventBus.getInstance().post(new DataLoadedEvent(frame));

    }

    /**
     * Remove out-of-view intervals from the cache.  This is O(N) where N = #frames X #intervals.   It is assumed
     * that N is small
     */
    private void trimCache() {
        List<AlignmentInterval> trimmedIntervals = new ArrayList<>();
        for(AlignmentInterval interval: intervalCache) {
            if (intervalInView(interval)) {
                trimmedIntervals.add(interval);
            }
        }
        intervalCache = trimmedIntervals;
    }


    private boolean intervalInView(AlignmentInterval interval) {

        for (ReferenceFrame frame : FrameManager.getFrames()) {
            if (interval.contains(frame.getCurrentRange())) {
                return true;
            }
        }
        return false;
    }


    AlignmentInterval loadInterval(String chr, int start, int end, AlignmentTrack.RenderOptions renderOptions) {

        String sequence = chrMappings.containsKey(chr) ? chrMappings.get(chr) : chr;

        DownsampleOptions downsampleOptions = new DownsampleOptions();

        final AlignmentTrack.BisulfiteContext bisulfiteContext =
                renderOptions != null ? renderOptions.bisulfiteContext : null;

        SpliceJunctionHelper spliceJunctionHelper = new SpliceJunctionHelper(this.loadOptions);

        ReadStats readStats = new ReadStats();

        AlignmentTileLoader.AlignmentTile t = checkReader().loadTile(sequence, start, end, spliceJunctionHelper,
                downsampleOptions, readStats, peStats, bisulfiteContext);
//
        if (getExperimentType() == null) {
            readStats.compute();
            inferType(readStats);
        }

        List<Alignment> alignments = t.getAlignments();
        List<DownsampledInterval> downsampledIntervals = t.getDownsampledIntervals();
        return new AlignmentInterval(chr, start, end, alignments, t.getCounts(), spliceJunctionHelper, downsampledIntervals);
    }

    /**
     * Some empirical metrics for determining experiment type
     *
     * @param readStats
     */
    private void inferType(ReadStats readStats) {

        if (readStats.readLengthStdDev > 100 || readStats.medianReadLength > 1000) {
            setExperimentType(AlignmentTrack.ExperimentType.THIRD_GEN);  // Could also use fracReadsWithIndels
        } else if (readStats.medianRefToReadRatio > 10) {
            setExperimentType(AlignmentTrack.ExperimentType.RNA);
        } else {
            setExperimentType(AlignmentTrack.ExperimentType.OTHER);
        }
    }

    private AlignmentTrack.ExperimentType getExperimentType() {
        return alignmentTrack == null ? null : alignmentTrack.getExperimentType();
    }

    private void setExperimentType(AlignmentTrack.ExperimentType type) {
        if (alignmentTrack != null) alignmentTrack.setExperimentType(type);
    }


    public  PackedAlignments getGroups(RenderContext context, AlignmentTrack.RenderOptions renderOptions) {
        AlignmentInterval interval = getLoadedInterval(context.getReferenceFrame());
        if (interval != null) {
            return interval.getPackedAlignments();
        } else {
            return null;
        }
    }

    public void clear() {
        intervalCache.clear();
    }

    public void dumpAlignments() {
        for (AlignmentInterval interval : intervalCache) {
            interval.dumpAlignments();
        }
    }

    /**
     * Find the first loaded interval for the specified chromosome and genomic {@code positon},
     * return the grouped alignments
     *
     * @param position
     * @param frame
     * @return alignmentRows, grouped and ordered by key
     */
    public PackedAlignments getGroupedAlignmentsContaining(double position, ReferenceFrame frame) {
        String chr = frame.getChrName();
        int start = (int) position;
        int end = start + 1;

        AlignmentInterval interval = getLoadedInterval(frame);
        if (interval == null) {
            return null;
        } else {
            PackedAlignments packedAlignments = interval.getPackedAlignments();
            if (packedAlignments != null && packedAlignments.contains(chr, start, end)) {
                return packedAlignments;
            } else {
                return null;
            }
        }
    }

    public int getNLevels() {
        int nLevels = 0;

        for (AlignmentInterval interval : intervalCache) {
            PackedAlignments packedAlignments = interval.getPackedAlignments();
            if (packedAlignments != null) {
                int intervalNLevels = packedAlignments.getNLevels();
                nLevels = Math.max(nLevels, intervalNLevels);
            }
        }
        return nLevels;
    }

    /**
     * Get the maximum group count among all the loaded intervals.  Normally there is one interval, but there
     * can be multiple if viewing split screen.
     */
    public int getMaxGroupCount() {
        int groupCount = 0;

        for (AlignmentInterval interval : intervalCache) {
            if (interval != null) {  // Not sure how this happens but it does
                PackedAlignments packedAlignments = interval.getPackedAlignments();
                if (packedAlignments != null) {
                    groupCount = Math.max(groupCount, packedAlignments.size());
                }
            }
        }
        return groupCount;
    }

    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        if (reader != null) {
            try {
                reader.close();
            } catch (IOException ex) {
                log.error("Error closing AlignmentQueryReader. ", ex);
            }
        }

    }

    public void updatePEStats(AlignmentTrack.RenderOptions renderOptions) {
        if (this.peStats != null) {
            for (PEStats stats : peStats.values()) {
                stats.computeInsertSize(renderOptions.getMinInsertSizePercentile(), renderOptions.getMaxInsertSizePercentile());
            }
        }
    }

    public SpliceJunctionHelper.LoadOptions getSpliceJunctionLoadOptions() {
        return loadOptions;
    }

    public void setMinJunctionCoverage(int minJunctionCoverage) {
        this.loadOptions = new SpliceJunctionHelper.LoadOptions(minJunctionCoverage, this.loadOptions.minReadFlankingWidth);
        for (AlignmentInterval interval : intervalCache) {
            interval.getSpliceJunctionHelper().setLoadOptions(this.loadOptions);
        }
    }

    public void alleleThresholdChanged() {
        coverageTrack.setSnpThreshold(PreferencesManager.getPreferences().getAsFloat(SAM_ALLELE_THRESHOLD));
    }

    public boolean isTenX() {
        return checkReader().isTenX();
    }

    public boolean isPhased() {
        return checkReader().isPhased();
    }

    public boolean isMoleculo() {
        return checkReader().isMoleculo();
    }

    public Collection<AlignmentInterval> getLoadedIntervals() {
        return intervalCache;
    }

    private AlignmentTileLoader checkReader() {
        try {
            String aPath = locator.getPath();
            if (AmazonUtils.isAwsS3Path(aPath) && !AmazonUtils.isS3PresignedValid(aPath)) {
                reader = new AlignmentTileLoader(AlignmentReaderFactory.getReader(locator));
            }
        } catch (MalformedURLException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return reader;
    }

    public static class DownsampleOptions {
        private boolean downsample;
        private int sampleWindowSize;
        private int maxReadCount;

        public DownsampleOptions() {
            IGVPreferences prefs = PreferencesManager.getPreferences();
            init(prefs.getAsBoolean(SAM_DOWNSAMPLE_READS),
                    prefs.getAsInt(SAM_SAMPLING_WINDOW),
                    prefs.getAsInt(SAM_SAMPLING_COUNT));
        }

        DownsampleOptions(boolean downsample, int sampleWindowSize, int maxReadCount) {
            init(downsample, sampleWindowSize, maxReadCount);
        }

        private void init(boolean downsample, int sampleWindowSize, int maxReadCount) {
            this.downsample = downsample;
            this.sampleWindowSize = sampleWindowSize;
            this.maxReadCount = maxReadCount;
        }

        public boolean isDownsample() {
            return downsample;
        }

        public int getSampleWindowSize() {
            return sampleWindowSize;
        }

        public int getMaxReadCount() {
            return maxReadCount;
        }

    }

    static class IntervalCache {

        private int maxSize;
        ArrayList<AlignmentInterval> intervals;

        public IntervalCache() {
            this(1);
        }

        public IntervalCache(int ms) {
            this.maxSize = Math.max(1, ms);
            intervals = new ArrayList<>(maxSize);
        }

        void setMaxSize(int ms, List<ReferenceFrame> frames) {
            this.maxSize = Math.max(1, ms);
            if (intervals.size() > maxSize) {
                // Reduce size.  Try to keep intervals that cover frame ranges.  This involves a linear search
                // of potentially (intervals.size X frames.size) elements.  Don't attempt if this number is too large
                if (frames.size() * intervals.size() < 25) {
                    ArrayList<AlignmentInterval> tmp = new ArrayList<>(maxSize);
                    for (AlignmentInterval interval : intervals) {
                        if (tmp.size() == maxSize) break;
                        for (ReferenceFrame frame : frames) {
                            Range range = frame.getCurrentRange();
                            if (interval.contains(range.getChr(), range.getStart(), range.getEnd())) {
                                tmp.add(interval);
                                break;
                            }
                        }
                    }
                    intervals = tmp;
                } else {
                    intervals = new ArrayList(intervals.subList(0, maxSize));
                    intervals.trimToSize();
                }
            }
        }

        public void add(AlignmentInterval interval) {
            if (intervals.size() >= maxSize) {
                intervals.remove(0);
            }
            intervals.add(interval);
        }

        public AlignmentInterval getIntervalForRange(Range range) {

            for (AlignmentInterval interval : intervals) {
                if (interval.contains(range.getChr(), range.getStart(), range.getEnd())) {
                    return interval;
                }
            }

            return null;

        }

        public Collection<AlignmentInterval> values() {
            return intervals;
        }

        public void clear() {
            intervals.clear();
        }
    }
}

