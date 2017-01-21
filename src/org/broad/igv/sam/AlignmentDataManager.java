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
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentTrack.SortOption;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.event.DataLoadedEvent;
import org.broad.igv.ui.event.IGVEventBus;
import org.broad.igv.ui.event.IGVEventObserver;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

public class AlignmentDataManager implements IGVEventObserver {

    private static Logger log = Logger.getLogger(AlignmentDataManager.class);

    private Map<ReferenceFrame, AlignmentInterval> intervalCache;

    private ResourceLocator locator;
    private HashMap<String, String> chrMappings = new HashMap();
    private volatile boolean isLoading = false;
    private AlignmentTileLoader reader;
    private CoverageTrack coverageTrack;
    private Map<String, PEStats> peStats;
    private ReadStats readStats;
    private AlignmentTrack.ExperimentType experimentType;
    private SpliceJunctionHelper.LoadOptions loadOptions;
    private Object loadLock = new Object();
    private boolean showAlignments = true;

    public AlignmentDataManager(ResourceLocator locator, Genome genome) throws IOException {
        this.locator = locator;
        reader = new AlignmentTileLoader(AlignmentReaderFactory.getReader(locator));
        readStats = new ReadStats(1000);
        peStats = new HashMap();
        initLoadOptions();
        initChrMap(genome);
        intervalCache = Collections.synchronizedMap(new HashMap<>());
        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
    }

    public void receiveEvent(Object event) {

        if (event instanceof FrameManager.ChangeEvent) {

            Collection<ReferenceFrame> frames = ((FrameManager.ChangeEvent) event).getFrames();
            Map<ReferenceFrame, AlignmentInterval> newCache = Collections.synchronizedMap(new HashMap<>());
            // Trim cache to include only current frames
            for(ReferenceFrame f : frames) {
                if(intervalCache.containsKey(f)) {
                    newCache.put(f, intervalCache.get(f));
                }
            }
            intervalCache = newCache;


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

    public void setExperimentType(AlignmentTrack.ExperimentType experimentType) {
        this.experimentType = experimentType;
    }

    public AlignmentTrack.ExperimentType getExperimentType() {
        return experimentType;
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

    public boolean hasIndex() {
        return reader.hasIndex();
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
        return intervalCache.get(frame);
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
        for (AlignmentInterval interval : intervalCache.values()) {
            interval.packAlignments(renderOptions);
        }
    }


    public boolean isLoaded(ReferenceFrame frame) {

        AlignmentInterval interval = intervalCache.get(frame);
        if (interval == null) {
            return false;
        } else {

            Range range = frame.getCurrentRange();

            return interval.contains(range.getChr(), range.getStart(), range.getEnd());
        }
    }


    public void load(ReferenceFrame referenceFrame,
                     AlignmentTrack.RenderOptions renderOptions,
                     boolean expandEnds) {

        if (isLoaded(referenceFrame)) return;  // Already loaded
        synchronized (loadLock) {
            final String chr = referenceFrame.getChrName();
            Range range = referenceFrame.getCurrentRange();

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
            loadAlignments(chr, adjustedStart, adjustedEnd, renderOptions, referenceFrame);
        }
    }

    private void loadAlignments(final String chr, final int start, final int end,
                                final AlignmentTrack.RenderOptions renderOptions,
                                final ReferenceFrame frame) {

        log.debug("Loading alignments: " + chr + ":" + start + "-" + end + " for " + AlignmentDataManager.this);

        AlignmentInterval loadedInterval = loadInterval(chr, start, end, renderOptions);
        intervalCache.put(frame, loadedInterval);

        packAlignments(renderOptions);
        IGVEventBus.getInstance().post(new DataLoadedEvent(frame));

        isLoading = false;


    }

    AlignmentInterval loadInterval(String chr, int start, int end, AlignmentTrack.RenderOptions renderOptions) {

        String sequence = chrMappings.containsKey(chr) ? chrMappings.get(chr) : chr;

        DownsampleOptions downsampleOptions = new DownsampleOptions();

        final AlignmentTrack.BisulfiteContext bisulfiteContext =
                renderOptions != null ? renderOptions.bisulfiteContext : null;

        SpliceJunctionHelper spliceJunctionHelper = new SpliceJunctionHelper(this.loadOptions);

        AlignmentTileLoader.AlignmentTile t = reader.loadTile(sequence, start, end, spliceJunctionHelper,
                downsampleOptions, readStats, peStats, bisulfiteContext, showAlignments);

        List<Alignment> alignments = t.getAlignments();
        List<DownsampledInterval> downsampledIntervals = t.getDownsampledIntervals();
        return new AlignmentInterval(chr, start, end, alignments, t.getCounts(), spliceJunctionHelper, downsampledIntervals);
    }


    public synchronized PackedAlignments getGroups(RenderContext context, AlignmentTrack.RenderOptions renderOptions) {
     //   load(context.getReferenceFrame(), renderOptions, false);
     //   Range range = context.getReferenceFrame().getCurrentRange();

        AlignmentInterval interval = intervalCache.get(context.getReferenceFrame());
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
        for (AlignmentInterval interval : intervalCache.values()) {
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

        AlignmentInterval interval = intervalCache.get(referenceFrame);
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

        for (AlignmentInterval interval : intervalCache.values()) {
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
        for (AlignmentInterval interval : intervalCache.values()) {
            if(interval != null) {  // Not sure how this happens but it does
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
        for (AlignmentInterval interval : intervalCache.values()) {
            interval.getSpliceJunctionHelper().setLoadOptions(this.loadOptions);
        }
    }

    public void alleleThresholdChanged() {
        coverageTrack.setSnpThreshold(PreferencesManager.getPreferences().getAsFloat(SAM_ALLELE_THRESHOLD));
    }

    public void setShowAlignments(boolean showAlignments) {
        this.showAlignments = showAlignments;
        clear();
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
        return intervalCache.values();
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

