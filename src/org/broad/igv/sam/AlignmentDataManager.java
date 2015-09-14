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

import com.google.common.eventbus.EventBus;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.AlignmentTrack.SortOption;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.event.DataLoadedEvent;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.*;

public class AlignmentDataManager implements IAlignmentDataManager {

    private static Logger log = Logger.getLogger(AlignmentDataManager.class);

    /**
     * Caches for loaded alignments and the relevant packing
     */
    private PositionCache<AlignmentInterval> loadedIntervalCache = new PositionCache<AlignmentInterval>();
    private PositionCache<PackedAlignments> packedAlignmentsCache = new PositionCache<PackedAlignments>();

    private HashMap<String, String> chrMappings = new HashMap();
    private volatile boolean isLoading = false;
    private AlignmentTileLoader reader;
    private CoverageTrack coverageTrack;

    private Map<String, PEStats> peStats;

    private AlignmentTrack.ExperimentType experimentType;

    private SpliceJunctionHelper.LoadOptions loadOptions;

    private Object loadLock = new Object();

    /**
     * This {@code EventBus} is typically used to notify listeners when new data
     * is loaded
     */
    private EventBus eventBus = new EventBus();

    private ResourceLocator locator;

    public AlignmentDataManager(ResourceLocator locator, Genome genome) throws IOException {
        this.locator = locator;
        reader = new AlignmentTileLoader(AlignmentReaderFactory.getReader(locator));
        peStats = new HashMap();
        initLoadOptions();
        initChrMap(genome);
    }

    void initLoadOptions() {
        this.loadOptions = new SpliceJunctionHelper.LoadOptions();
    }

    /**
     * Create an alias -> chromosome lookup map.  Enable loading BAM files that use alternative names for chromosomes,
     * provided the alias has been defined  (e.g. 1 -> chr1,  etc).
     */
    private void initChrMap(Genome genome) {
        if (genome != null) {
            List<String> seqNames = reader.getSequenceNames();
            if (seqNames != null) {
                for (String chr : seqNames) {
                    String alias = genome.getChromosomeAlias(chr);
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

    /**
     * The set of sequences found in the file.
     * May be null
     *
     * @return
     */
    public List<String> getSequenceNames() {
        return reader.getSequenceNames();
    }


    public boolean isIonTorrent() {
        Set<String> platforms = reader.getPlatforms();
        if (platforms != null) {
            return platforms.contains("IONTORRENT");
        }
        return false;
    }

    public Collection<AlignmentInterval> getLoadedIntervals() {
        return this.loadedIntervalCache.values();
    }

    public AlignmentInterval getLoadedInterval(Range range) {
        return loadedIntervalCache.getForRange(range);
    }

    /**
     * Sort rows group by group
     *
     * @param option
     * @param location
     */
    public boolean sortRows(SortOption option, ReferenceFrame frame, double location, String tag) {
        PackedAlignments packedAlignments = packedAlignmentsCache.getForRange(frame.getCurrentRange());
        AlignmentInterval interval = loadedIntervalCache.getForRange(frame.getCurrentRange());
        if (packedAlignments == null || interval == null) {
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
    boolean packAlignments(AlignmentTrack.RenderOptions renderOptions) {

        List<ReferenceFrame> frameList = FrameManager.getFrames();
        List<AlignmentInterval> intervalList = new ArrayList<AlignmentInterval>(frameList.size());

        this.packedAlignmentsCache.clear();
        this.packedAlignmentsCache.setMaxEntries(2 * intervalList.size());

        for (ReferenceFrame frame : frameList) {
            AlignmentInterval interval = loadedIntervalCache.getForRange(frame.getCurrentRange());

            if (interval == null) {
                return false;
            }

            final AlignmentPacker alignmentPacker = new AlignmentPacker();
            PackedAlignments packedAlignments = alignmentPacker.packAlignments(interval, renderOptions);

            //We cache by the interval range because this will generally be buffered/expanded, whereas the frame
            //will be to-the-pixel (meaning a slight scroll triggers a repack

            this.packedAlignmentsCache.put(interval.getRange(), packedAlignments);
        }

        return true;
    }

    public void load(RenderContext context,
                     AlignmentTrack.RenderOptions renderOptions,
                     boolean expandEnds) {

        synchronized (loadLock) {
            final String chr = context.getChr();
            final int start = (int) context.getOrigin();
            final int end = (int) context.getEndLocation();
            AlignmentInterval loadedInterval = loadedIntervalCache.getForRange(context.getReferenceFrame().getCurrentRange());

            int adjustedStart = start;
            int adjustedEnd = end;
            // Expand the interval by the lesser of  +/- a 2 screens, or max visible range
            int windowSize = Math.min(4*(end - start), PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_MAX_VISIBLE_RANGE) * 1000);
            int center = (end + start) / 2;
            int expand = Math.max(end - start, windowSize / 2);

            if (loadedInterval != null) {
                // First see if we have a loaded interval that fully contain the requested interval.
                // If so, we don't need to load it
                if (loadedInterval.contains(chr, start, end)) {
                    return;
                }
            }

            if (expandEnds) {
                adjustedStart = Math.max(0, Math.min(start, center - expand));
                adjustedEnd = Math.max(end, center + expand);
            }
            loadAlignments(chr, adjustedStart, adjustedEnd, renderOptions, context);
        }

    }

    public synchronized PackedAlignments getGroups(RenderContext context, AlignmentTrack.RenderOptions renderOptions) {
        load(context, renderOptions, false);
        Range range = context.getReferenceFrame().getCurrentRange();
        if (!packedAlignmentsCache.containsRange(range)) {
            packAlignments(renderOptions);
        }
        return packedAlignmentsCache.getForRange(context.getReferenceFrame().getCurrentRange());
    }

    public void clear() {
        // reader.clearCache();
        loadedIntervalCache.clear();
        packedAlignmentsCache.clear();
    }

    public synchronized void loadAlignments(final String chr, final int start, final int end,
                                            final AlignmentTrack.RenderOptions renderOptions,
                                            final RenderContext context) {

        if (isLoading || chr.equals(Globals.CHR_ALL)) {
            return;
        }
        loadedIntervalCache.setMaxEntries(2 * FrameManager.getFrames().size());
        isLoading = true;

        NamedRunnable runnable = new NamedRunnable() {

            public String getName() {
                return "loadAlignments";
            }

            public void run() {

                log.debug("Loading alignments: " + chr + ":" + start + "-" + end + " for " + AlignmentDataManager.this);

                AlignmentInterval loadedInterval = loadInterval(chr, start, end, renderOptions);
                loadedIntervalCache.put(loadedInterval.getRange(), loadedInterval);

                List<ReferenceFrame> frameList = context != null ? Arrays.asList(context.getReferenceFrame()) : null;
                packAlignments(renderOptions);
                getEventBus().post(new DataLoadedEvent(context));

                isLoading = false;
            }
        };
        LongRunningTask.submit(runnable);
    }

    AlignmentInterval loadInterval(String chr, int start, int end, AlignmentTrack.RenderOptions renderOptions) {

        String sequence = chrMappings.containsKey(chr) ? chrMappings.get(chr) : chr;

        DownsampleOptions downsampleOptions = new DownsampleOptions();

        final AlignmentTrack.BisulfiteContext bisulfiteContext =
                renderOptions != null ? renderOptions.bisulfiteContext : null;

        ProgressMonitor monitor = null;
        //Show cancel button
        if (IGV.hasInstance() && !Globals.isBatch() && !Globals.isHeadless()) {
            ActionListener cancelListener = new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    AlignmentTileLoader.cancelReaders();
                }
            };
            IGV.getInstance().getContentPane().getStatusBar().activateCancelButton(cancelListener);
        }

        SpliceJunctionHelper spliceJunctionHelper = new SpliceJunctionHelper(this.loadOptions);
        AlignmentTileLoader.AlignmentTile t = reader.loadTile(sequence, start, end, spliceJunctionHelper,
                downsampleOptions, peStats, bisulfiteContext, monitor);

        List<Alignment> alignments = t.getAlignments();
        List<DownsampledInterval> downsampledIntervals = t.getDownsampledIntervals();
        return new AlignmentInterval(chr, start, end, alignments, t.getCounts(), spliceJunctionHelper, downsampledIntervals);
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

        PackedAlignments packedAlignments = packedAlignmentsCache.getForRange(referenceFrame.getCurrentRange());
        if (packedAlignments != null && packedAlignments.contains(chr, start, end)) {
            return packedAlignments;
        }
        return null;
    }

    public int getNLevels() {
        int nLevels = 0;
        for (PackedAlignments packedAlignments : packedAlignmentsCache.values()) {
            int intervalNLevels = packedAlignments.getNLevels();
            nLevels = Math.max(nLevels, intervalNLevels);
        }
        return nLevels;
    }

    /**
     * Get the maximum group count among all the loaded intervals.  Normally there is one interval, but there
     * can be multiple if viewing split screen.
     */
    public int getMaxGroupCount() {
        int groupCount = 0;
        for (PackedAlignments packedAlignments : packedAlignmentsCache.values()) {
            groupCount = Math.max(groupCount, packedAlignments.size());
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
                stats.compute(renderOptions.getMinInsertSizePercentile(), renderOptions.getMaxInsertSizePercentile());
            }
        }
    }

    public EventBus getEventBus() {
        return eventBus;
    }

    public SpliceJunctionHelper.LoadOptions getSpliceJunctionLoadOptions() {
        return loadOptions;
    }

    public void setMinJunctionCoverage(int minJunctionCoverage) {
        this.loadOptions = new SpliceJunctionHelper.LoadOptions(minJunctionCoverage, this.loadOptions.minReadFlankingWidth);
        for (AlignmentInterval interval : getLoadedIntervals()) {
            interval.getSpliceJunctionHelper().setLoadOptions(this.loadOptions);
        }
    }

    PositionCache getCache() {
        return this.loadedIntervalCache;
    }

    public void alleleThresholdChanged() {
        coverageTrack.setSnpThreshold(PreferenceManager.getInstance().getAsFloat(PreferenceManager.SAM_ALLELE_THRESHOLD));
    }

    public static class DownsampleOptions {
        private boolean downsample;
        private int sampleWindowSize;
        private int maxReadCount;

        public DownsampleOptions() {
            PreferenceManager prefs = PreferenceManager.getInstance();
            init(prefs.getAsBoolean(PreferenceManager.SAM_DOWNSAMPLE_READS),
                    prefs.getAsInt(PreferenceManager.SAM_SAMPLING_WINDOW),
                    prefs.getAsInt(PreferenceManager.SAM_SAMPLING_COUNT));
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
}

