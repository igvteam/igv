/*
 * Copyright (c) 2007-2014 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
package org.broad.igv.sam;

import com.google.common.eventbus.EventBus;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.AlignmentTrack.SortOption;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.event.DataLoadedEvent;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.ArrayHeapObjectSorter;
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
     * Map of reference frame name -> alignment interval
     */
    private Map<String, AlignmentInterval> loadedIntervalMap = new HashMap<String, AlignmentInterval>();

    private Map<String, PackedAlignments> packedAlignmentsMap = new HashMap<String, PackedAlignments>();

    private HashMap<String, String> chrMappings = new HashMap();
    private volatile boolean isLoading = false;
    private AlignmentTileLoader reader;
    private CoverageTrack coverageTrack;

    private static final int MAX_ROWS = 1000000;
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

    public Collection<AlignmentInterval> getAllLoadedIntervals() {
        return loadedIntervalMap.values();
    }

    Collection<String> getLoadedIntervalNames() {
        return loadedIntervalMap.keySet();
    }

    /**
     * Return the loaded interval for the specified frame (by name).  Note this can be null if the interval isn't loaded
     * yet.
     *
     * @param frameName
     * @return
     */
    public AlignmentInterval getLoadedInterval(String frameName) {
        return loadedIntervalMap.get(frameName);
    }

    /**
     * Sort rows group by group
     *
     * @param option
     * @param location
     */
    public void sortRows(SortOption option, String frameName, double location, String tag) {
        PackedAlignments packedAlignments = packedAlignmentsMap.get(frameName);
        AlignmentInterval interval = loadedIntervalMap.get(frameName);
        if (packedAlignments == null || interval == null) {
            return;
        }

        for (List<Row> alignmentRows : packedAlignments.values()) {
            for (Row row : alignmentRows) {
                row.updateScore(option, location, interval, tag);
            }
            Collections.sort(alignmentRows);
        }
    }

    public void setViewAsPairs(boolean option, AlignmentTrack.RenderOptions renderOptions) {
        if (option == renderOptions.isViewPairs()) {
            return;
        }

        boolean currentPairState = renderOptions.isViewPairs();
        renderOptions.setViewPairs(option);

        for (ReferenceFrame frame : FrameManager.getFrames()) {
            repackAlignments(frame.getName(), currentPairState, renderOptions);
        }
    }

    private void repackAlignments(String frameName, boolean currentPairState, AlignmentTrack.RenderOptions renderOptions) {

        if (currentPairState) {
            PackedAlignments packedAlignments = packedAlignmentsMap.get(frameName);
            if (packedAlignments == null) {
                return;
            }

            List<Alignment> alignments = new ArrayList<Alignment>(Math.min(50000, packedAlignments.size() * 10000));
            int intervalEnd = -1;
            for (List<Row> alignmentRows : packedAlignments.values()) {
                for (Row row : alignmentRows) {
                    for (Alignment al : row.alignments) {
                        intervalEnd = Math.max(intervalEnd, al.getEnd());
                        if (al instanceof PairedAlignment) {
                            PairedAlignment pair = (PairedAlignment) al;
                            alignments.add(pair.firstAlignment);
                            if (pair.secondAlignment != null) {
                                alignments.add(pair.secondAlignment);
                            }
                        } else {
                            alignments.add(al);
                        }
                    }
                }
            }


            // ArrayHeapObjectSorter sorts in place (no additional memory required).
            ArrayHeapObjectSorter<Alignment> heapSorter = new ArrayHeapObjectSorter();
            heapSorter.sort(alignments, new Comparator<Alignment>() {
                public int compare(Alignment alignment, Alignment alignment1) {
                    return alignment.getStart() - alignment1.getStart();
                }
            });

            PackedAlignments tmp = (new AlignmentPacker()).packAlignments(
                    alignments.iterator(),
                    intervalEnd,
                    renderOptions);

            packedAlignmentsMap.put(frameName, tmp);

        } else {
            repackAlignments(frameName, renderOptions);
        }
    }

    /**
     * Repack currently loaded alignments of the provided reference frame
     *
     * @param frameName
     * @param renderOptions
     * @see AlignmentPacker#packAlignments(java.util.Iterator, int, org.broad.igv.sam.AlignmentTrack.RenderOptions)
     */
    public void repackAlignments(String frameName, AlignmentTrack.RenderOptions renderOptions) {

        AlignmentInterval loadedInterval = loadedIntervalMap.get(frameName);
        if (loadedInterval == null) {
            return;
        }

        Iterator<Alignment> iter = loadedInterval.getAlignmentIterator();
        PackedAlignments packedAlignments = (new AlignmentPacker()).packAlignments(
                iter,
                loadedInterval.getEnd(),
                renderOptions);

        this.packedAlignmentsMap.put(frameName, packedAlignments);
    }


    public void load(RenderContext context,
                     AlignmentTrack.RenderOptions renderOptions,
                     boolean expandEnds) {

        synchronized (loadLock) {
            final String chr = context.getChr();
            final int start = (int) context.getOrigin();
            final int end = (int) context.getEndLocation();
            AlignmentInterval loadedInterval = loadedIntervalMap.get(context.getReferenceFrame().getName());

            int adjustedStart = start;
            int adjustedEnd = end;
            int windowSize = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_MAX_VISIBLE_RANGE) * 1000;
            int center = (end + start) / 2;
            int expand = Math.max(end - start, windowSize / 2);

            if (loadedInterval != null) {
                // First see if we have a loaded interval that fully contain the requested interval.  If yes we're done
                if (loadedInterval.contains(chr, start, end)) {
                    // Requested interval is fully contained in the existing one, we're done
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

    public synchronized PackedAlignments getGroups(RenderContext context,
                                                                                     AlignmentTrack.RenderOptions renderOptions) {

        load(context, renderOptions, false);
        return packedAlignmentsMap.get(context.getReferenceFrame().getName());
    }

    public void clear() {
        // reader.clearCache();
        loadedIntervalMap.clear();
    }

    public synchronized void loadAlignments(final String chr, final int start, final int end,
                                            final AlignmentTrack.RenderOptions renderOptions,
                                            final RenderContext context) {

        if (isLoading || chr.equals(Globals.CHR_ALL)) {
            return;
        }

        isLoading = true;

        NamedRunnable runnable = new NamedRunnable() {

            public String getName() {
                return "loadAlignments";
            }

            public void run() {

                log.debug("Loading alignments: " + chr + ":" + start + "-" + end + " for " + AlignmentDataManager.this);

                AlignmentInterval loadedInterval = loadInterval(chr, start, end, renderOptions);
                final AlignmentPacker alignmentPacker = new AlignmentPacker();
                PackedAlignments packedAlignments = alignmentPacker.packAlignments(loadedInterval.getAlignmentIterator(), end, renderOptions);

                ReferenceFrame frame = context != null ? context.getReferenceFrame() : null;
                addLoadedInterval(frame, loadedInterval, packedAlignments);

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
        return new AlignmentInterval(chr, start, end, alignments, t.getCounts(), spliceJunctionHelper, downsampledIntervals, renderOptions);
    }

    private void addLoadedInterval(ReferenceFrame frame, AlignmentInterval interval, PackedAlignments packedAlignments) {
        String frameName = frame != null ? frame.getName() : FrameManager.DEFAULT_FRAME_NAME;
        loadedIntervalMap.put(frameName, interval);
        packedAlignmentsMap.put(frameName, packedAlignments);
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
        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        if (loadedInterval == null) return null;

        PackedAlignments packedAlignments = packedAlignmentsMap.get(referenceFrame.getName());
        if (packedAlignments != null && loadedInterval.contains(chr, start, end)) {
            return packedAlignments;
        }
        return null;
    }

    public int getNLevels() {
        int nLevels = 0;
        for (PackedAlignments packedAlignments : packedAlignmentsMap.values()) {
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
        for (PackedAlignments packedAlignments : packedAlignmentsMap.values()) {
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
        for (AlignmentInterval interval : getAllLoadedIntervals()) {
            interval.getSpliceJunctionHelper().setLoadOptions(this.loadOptions);
        }
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

