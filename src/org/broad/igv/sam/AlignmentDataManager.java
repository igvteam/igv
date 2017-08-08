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
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentTrack.SortOption;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.track.RenderContext;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * Manages data loading for a single alignment file.  Shared between alignment, coverage, and junction
 * tracks.
 */

public class AlignmentDataManager implements IGVEventObserver {

    private static Logger log = Logger.getLogger(AlignmentDataManager.class);


    private Collection<AlignmentInterval> intervalCache;
    private ResourceLocator locator;
    private HashMap<String, String> chrMappings = new HashMap();
    private Set<Range> isLoading = new HashSet<>();
    private AlignmentTileLoader reader;
    private CoverageTrack coverageTrack;
    private Map<String, PEStats> peStats;
    private SpliceJunctionHelper.LoadOptions loadOptions;
    private Object loadLock = new Object();
    private boolean showAlignments = true;
    private AlignmentTrack.ExperimentType inferredExperimentType;

    public AlignmentDataManager(ResourceLocator locator, Genome genome) throws IOException {
        this.locator = locator;
        reader = new AlignmentTileLoader(AlignmentReaderFactory.getReader(locator));
        peStats = new HashMap();
        initLoadOptions();
        initChrMap(genome);
        intervalCache = Collections.synchronizedList(new ArrayList<>());

        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(RefreshEvent.class, this);
    }

    public void receiveEvent(Object event) {

        if (event instanceof FrameManager.ChangeEvent) {

            Collection<ReferenceFrame> frames = ((FrameManager.ChangeEvent) event).getFrames();
            Collection<AlignmentInterval> newCache = Collections.synchronizedList(new ArrayList<>());

            // Trim cache to include only current frames

            for (ReferenceFrame f : frames) {
                AlignmentInterval i = getLoadedInterval(f);
                if (i != null) {
                    newCache.add(i);
                }
            }
            intervalCache = newCache;


        } else if (event instanceof RefreshEvent) {
            clear();
        } else {
            log.info("Unknown event type: " + event.getClass());
        }
    }

    void initLoadOptions() {
        this.loadOptions = new SpliceJunctionHelper.LoadOptions();
    }

    /**
     * Create an alias -> chromosome lookup map.  Enable loading BAM files that use alternative names for chromosomes,
     * provided the alias has been defined  (e.g. 1 -> chr1,  etc).
     */
    private void initChrMap(Genome genome) throws IOException {
        if (genome != null) {
            List<String> seqNames = reader.getSequenceNames();
            if (seqNames != null) {
                for (String chr : seqNames) {
                    String alias = genome.getCanonicalChrName(chr);
                    chrMappings.put(alias, chr);
                }
            }
        }
    }

    public AlignmentTileLoader getReader() {
        return reader;
    }

    public ResourceLocator getLocator() {
        return locator;
    }

    public Map<String, PEStats> getPEStats() {
        return peStats;
    }

    public boolean isPairedEnd() {
        return reader.isPairedEnd();
    }

    public boolean hasYCTags() {
        return reader.hasYCTags();
    }

    public boolean hasIndex() {
        return reader.hasIndex();
    }

    public void setInferredExperimentType(AlignmentTrack.ExperimentType inferredExperimentType) {
        if (inferredExperimentType != this.inferredExperimentType) {
            ExperimentTypeChangeEvent event = new ExperimentTypeChangeEvent(this, inferredExperimentType);
            this.inferredExperimentType = inferredExperimentType;
            IGVEventBus.getInstance().post(event);
        }
    }

    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public double getMinVisibleScale() {
        IGVPreferences prefs = PreferencesManager.getPreferences();
        float maxRange = prefs.getAsFloat(SAM_MAX_VISIBLE_RANGE);
        return (maxRange * 1000) / 700;
    }


    /**
     * The set of sequences found in the file.
     * May be null
     *
     * @return
     */
    public List<String> getSequenceNames() throws IOException {
        return reader.getSequenceNames();
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

    public boolean isLoading(ReferenceFrame frame) {

        Range range = frame.getCurrentRange();
        for (Range r : isLoading) {
            if (r.contains(range)) return true;
        }
        return false;
    }


    public void load(ReferenceFrame referenceFrame,
                     AlignmentTrack.RenderOptions renderOptions,
                     boolean expandEnds) {

        if (isLoaded(referenceFrame)) return;  // Already loaded

        if (isLoading(referenceFrame)) return;   // Already oading

        synchronized (loadLock) {
            Range range = referenceFrame.getCurrentRange();

            isLoading.add(range);

            final String chr = referenceFrame.getChrName();

            final int start = (int) range.getStart();
            final int end = (int) range.getEnd();
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


            log.debug("Loading alignments: " + chr + ":" + adjustedStart + "-" + adjustedEnd + " for " + AlignmentDataManager.this);

            AlignmentInterval loadedInterval = loadInterval(chr, adjustedStart, adjustedEnd, renderOptions);
            intervalCache.add(loadedInterval);

            packAlignments(renderOptions);
            isLoading.remove(range);

            //  IGVEventBus.getInstance().post(new DataLoadedEvent(referenceFrame));

        }
    }


    AlignmentInterval loadInterval(String chr, int start, int end, AlignmentTrack.RenderOptions renderOptions) {

        String sequence = chrMappings.containsKey(chr) ? chrMappings.get(chr) : chr;

        DownsampleOptions downsampleOptions = new DownsampleOptions();

        final AlignmentTrack.BisulfiteContext bisulfiteContext =
                renderOptions != null ? renderOptions.bisulfiteContext : null;

        SpliceJunctionHelper spliceJunctionHelper = new SpliceJunctionHelper(this.loadOptions);

        ReadStats readStats = new ReadStats();

        AlignmentTileLoader.AlignmentTile t = reader.loadTile(sequence, start, end, spliceJunctionHelper,
                downsampleOptions, readStats, peStats, bisulfiteContext, showAlignments);

        if (inferredExperimentType == null && !Globals.VERSION.contains("2.4")) {
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
            setInferredExperimentType(AlignmentTrack.ExperimentType.THIRD_GEN);  // Could also use fracReadsWithIndels
        } else if (readStats.medianRefToReadRatio > 10) {
            setInferredExperimentType(AlignmentTrack.ExperimentType.RNA);
        } else {
            setInferredExperimentType(AlignmentTrack.ExperimentType.OTHER);
        }
    }


    public synchronized PackedAlignments getGroups(RenderContext context, AlignmentTrack.RenderOptions renderOptions) {
        //   load(context.getReferenceFrame(), renderOptions, false);
        //   Range range = context.getReferenceFrame().getCurrentRange();

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
     * @param referenceFrame
     * @return alignmentRows, grouped and ordered by key
     */
    public PackedAlignments getGroupedAlignmentsContaining(double position, ReferenceFrame referenceFrame) {
        String chr = referenceFrame.getChrName();
        int start = (int) position;
        int end = start + 1;

        AlignmentInterval interval = getLoadedInterval(referenceFrame);
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

    public void setShowAlignments(boolean showAlignments) {
        if (showAlignments != this.showAlignments) {
            this.showAlignments = showAlignments;
            if (showAlignments == false) {
                dumpAlignments();
            } else {
                // Change from false => true,  need to reload
                intervalCache.clear();
            }
        }

    }

    public boolean isTenX() {
        return reader.isTenX();
    }

    public boolean isPhased() {
        return reader.isPhased();
    }

    public boolean isMoleculo() {
        return reader.isMoleculo();
    }

    public Collection<AlignmentInterval> getLoadedIntervals() {
        return intervalCache;
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

