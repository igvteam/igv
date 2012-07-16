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
import org.broad.igv.data.Interval;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.AlignmentTrack.SortOption;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ArrayHeapObjectSorter;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.CachedIntervals;

import java.io.IOException;
import java.util.*;

public class AlignmentDataManager {

    private static Logger log = Logger.getLogger(AlignmentDataManager.class);

    private static final int DEFAULT_DEPTH = 10;

    /**
     * Map of reference frame -> alignment interval
     */
    //TODO -- this is a  potential memory leak, this map needs cleared when the gene list changes
    private CachedIntervals<AlignmentInterval> loadedIntervalMap = new CachedIntervals(CACHE_SIZE, (int) 1e6);

    private HashMap<String, String> chrMappings = new HashMap();
    private boolean isLoading = false;
    private CachingQueryReader reader;
    private CoverageTrack coverageTrack;

    private static final int MAX_ROWS = 1000000;
    private Map<String, PEStats> peStats;

    private AlignmentTrack.ExperimentType experimentType;

    private boolean showSpliceJunctions;

    static final int CACHE_SIZE = 5;

    /**
     * Don't wont to allow intervals to get too big. In general,
     * will trim a cached interval it's bigger than MAX_INTERVAL_MULTIPLE * current interval size
     */
    static final int MAX_INTERVAL_MULTIPLE = 5;


    public AlignmentDataManager(ResourceLocator locator, Genome genome) throws IOException {

        PreferenceManager prefs = PreferenceManager.getInstance();
        reader = new CachingQueryReader(AlignmentReaderFactory.getReader(locator));
        peStats = new HashMap();
        showSpliceJunctions = prefs.getAsBoolean(PreferenceManager.SAM_SHOW_JUNCTION_TRACK);
        initChrMap(genome);
    }

    public void setShowSpliceJunctions(boolean showSpliceJunctions) {
        this.showSpliceJunctions = showSpliceJunctions;
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
        if (experimentType == AlignmentTrack.ExperimentType.BISULFITE) {
            showSpliceJunctions = false;
        } else {
            showSpliceJunctions = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_JUNCTION_TRACK);
        }
    }

    public AlignmentTrack.ExperimentType getExperimentType() {
        return experimentType;
    }

    public CachingQueryReader getReader() {
        return reader;
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


    /**
     * Return the loaded interval for the specified frame.  Note this can be null if the interval isn't loaded
     * yet.
     *
     * @param frame
     * @return
     */
    public Collection<AlignmentInterval> getLoadedIntervals(ReferenceFrame frame) {
        return loadedIntervalMap.get(frame.getChrName());
    }


    public void sortRows(SortOption option, ReferenceFrame referenceFrame, double location, String tag) {
        List<AlignmentInterval> loadedIntervals = loadedIntervalMap.get(referenceFrame.getChrName());
        if(loadedIntervals == null) return;
        for (AlignmentInterval loadedInterval : loadedIntervals) {
            if (loadedInterval != null) {
                loadedInterval.sortRows(option, location, tag);
            }
        }
    }

    public void setViewAsPairs(boolean option, AlignmentTrack.RenderOptions renderOptions) {
        if (option == renderOptions.isViewPairs()) {
            return;
        }

        boolean currentPairState = renderOptions.isViewPairs();
        renderOptions.setViewPairs(option);

        for (ReferenceFrame frame : FrameManager.getFrames()) {
            repackAlignments(frame, currentPairState, renderOptions);
        }
    }

    private void repackAlignments(ReferenceFrame referenceFrame, boolean currentPairState, AlignmentTrack.RenderOptions renderOptions) {

        if (currentPairState == true) {
            List<AlignmentInterval> loadedIntervals = loadedIntervalMap.get(referenceFrame.getChrName());
            if (loadedIntervals == null) {
                return;
            }


            for (AlignmentInterval loadedInterval : loadedIntervals) {
                Map<String, List<AlignmentInterval.Row>> groupedAlignments = loadedInterval.getGroupedAlignments();
                List<Alignment> alignments = new ArrayList(Math.min(50000, groupedAlignments.size() * 10000));
                for (List<AlignmentInterval.Row> alignmentRows : groupedAlignments.values()) {
                    for (AlignmentInterval.Row row : alignmentRows) {
                        for (Alignment al : row.alignments) {
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

                // When repacking keep all currently loaded alignments (don't limit to levels)
                int max = Integer.MAX_VALUE;
                LinkedHashMap<String, List<AlignmentInterval.Row>> tmp = (new AlignmentPacker()).packAlignments(
                        alignments.iterator(),
                        loadedInterval.getEnd(),
                        renderOptions);

                loadedInterval.setAlignmentRows(tmp, renderOptions);
            }
        } else {
            repackAlignments(referenceFrame, renderOptions);
        }
    }

    /**
     * Repack currently loaded alignments.
     *
     * @param referenceFrame
     */
    public void repackAlignments(ReferenceFrame referenceFrame, AlignmentTrack.RenderOptions renderOptions) {

        List<AlignmentInterval> loadedIntervals = loadedIntervalMap.get(referenceFrame.getChrName());
        if (loadedIntervals == null) {
            return;
        }


        for (AlignmentInterval loadedInterval : loadedIntervals) {
            Iterator<Alignment> iter = loadedInterval.getAlignmentIterator();

            // When repacking keep all currently loaded alignments (don't limit to levels)
            int max = Integer.MAX_VALUE;
            LinkedHashMap<String, List<AlignmentInterval.Row>> alignmentRows = (new AlignmentPacker()).packAlignments(
                    iter,
                    loadedInterval.getEnd(),
                    renderOptions);

            loadedInterval.setAlignmentRows(alignmentRows, renderOptions);
        }
    }

    public synchronized void preload(RenderContext context,
                                     AlignmentTrack.RenderOptions renderOptions,
                                     AlignmentTrack.BisulfiteContext bisulfiteContext,
                                     boolean expandEnds) {

        final String chr = context.getChr();
        final int start = (int) context.getOrigin();
        final int end = (int) context.getEndLocation();
        List<AlignmentInterval> loadedIntervals = loadedIntervalMap.get(context.getReferenceFrame().getChrName());
        boolean haveInterval = false;

        int adjustedStart = start;
        int adjustedEnd = end;

        if (loadedIntervals != null) {
            for (AlignmentInterval loadedInterval : loadedIntervals) {
                if (loadedInterval.contains(chr, start, end)) {
                    haveInterval = true;
                    break;
                } else if (loadedInterval.overlaps(chr, start, end, context.getZoom())) {
                    //We have part of the interval, only need to lead the portion
                    //we don't have
                    adjustedStart = Math.max(start, loadedInterval.getStart() - 1);
                    adjustedEnd = Math.min(end, loadedInterval.getEnd() + 1);
                }
            }
        }


        if (expandEnds) {
            int length = Math.max(100000, end - start);
            adjustedStart = Math.max(0, start - length / 2);
            adjustedEnd = end + length / 2;
        }

        // If we've moved out of the loaded interval start a new load.
        if (!haveInterval) {
            loadAlignments(chr, adjustedStart, adjustedEnd, renderOptions, context, bisulfiteContext);
        }

    }

    public synchronized LinkedHashMap<String, List<AlignmentInterval.Row>> getGroups(RenderContext context,
                                                                                     AlignmentTrack.RenderOptions renderOptions,
                                                                                     AlignmentTrack.BisulfiteContext bisulfiteContext) {

        final String chr = context.getChr();
        final int start = (int) context.getOrigin();
        final int end = (int) context.getEndLocation();

        preload(context, renderOptions, bisulfiteContext, false);

        List<AlignmentInterval> overlaps = loadedIntervalMap.getOverlaps(chr, start, end, context.getZoom());
        if (overlaps != null && overlaps.size() >= 1){
            // If there is any overlap in the loaded interval and the requested interval return it.
            return overlaps.get(0).getGroupedAlignments();
        }
        return null;
    }


    public void clear() {
        reader.clearCache();
        loadedIntervalMap.clear();
    }

    public synchronized void loadAlignments(final String chr, final int start, final int end,
                                            final AlignmentTrack.RenderOptions renderOptions,
                                            final RenderContext context,
                                            final AlignmentTrack.BisulfiteContext bisulfiteContext) {

        if (isLoading || chr.equals(Globals.CHR_ALL)) {
            return;
        }

        log.debug("Load alignments.  isLoading=" + isLoading);
        isLoading = true;
        NamedRunnable runnable = new NamedRunnable() {

            public String getName() {
                return "loadAlignments";
            }

            public void run() {

                log.debug("Loading alignments: " + chr + ":" + start + "-" + end);

                // Expand start and end to facilitate panning

                int expandLength = reader.getTileSize(chr) / 2;
                int intervalStart = start - expandLength;
                int intervalEnd = end + expandLength;

                DownsampleOptions downsampleOptions = new DownsampleOptions();

                CloseableIterator<Alignment> iter = null;
                //try {

                String sequence = chrMappings.containsKey(chr) ? chrMappings.get(chr) : chr;

                List<AlignmentCounts> counts = new ArrayList();
                List<CachingQueryReader.DownsampledInterval> downsampledIntervals = new ArrayList<CachingQueryReader.DownsampledInterval>();
                List<SpliceJunctionFeature> spliceJunctions = null;
                if (showSpliceJunctions) {
                    spliceJunctions = new ArrayList<SpliceJunctionFeature>();
                }

                iter = reader.query(sequence, intervalStart, intervalEnd, counts, spliceJunctions, downsampledIntervals,
                        downsampleOptions, peStats, bisulfiteContext);

                final AlignmentPacker alignmentPacker = new AlignmentPacker();

                LinkedHashMap<String, List<AlignmentInterval.Row>> alignmentRows = alignmentPacker.packAlignments(iter,
                        intervalEnd, renderOptions);

                AlignmentInterval loadedInterval = new AlignmentInterval(chr, intervalStart, intervalEnd,
                        alignmentRows, counts, spliceJunctions, downsampledIntervals, renderOptions);

                addLoadedInterval(context, loadedInterval);


                if (coverageTrack != null) {
                    coverageTrack.rescale(context.getReferenceFrame());
                }

                // TODO --- we need to force a repaint of the coverageTrack, which might not be in the same panel
                if (context.getPanel() != null) {
                    context.getPanel().invalidate();
                    context.getPanel().repaint();
                }

//                } catch (Exception exception) {
//                    log.error("Error loading alignments", exception);
//                    JOptionPane.showMessageDialog(IGV.getMainFrame(), "Error reading file: " + exception.getMessage());
//                } finally {
//
                if (iter != null) {
                    iter.close();
                }
                isLoading = false;
//                }
            }
        };

        LongRunningTask.submit(runnable);


    }

    private void addLoadedInterval(RenderContext context, AlignmentInterval interval) {
        loadedIntervalMap.setMaxIntervalSize(MAX_INTERVAL_MULTIPLE * (int) (context.getEndLocation() - context.getOrigin()));
        loadedIntervalMap.put(interval);
//        String key = context.getReferenceFrame().getChrName();
//        List<AlignmentInterval> currentValue = loadedIntervalMap.get(key);
//        if (currentValue != null) {
//            for (AlignmentInterval loadedInterval : currentValue) {
//                if (loadedInterval.overlaps(interval.getChr(), interval.getStart(), interval.getEnd(), interval.getZoom())) {
//                    loadedInterval.merge(interval);
//
//                    //Prevent interval from growing without bound
//                    int intervalSize = loadedInterval.getEnd() - loadedInterval.getStart();
//                    if (intervalSize > MAX_INTERVAL_MULTIPLE * (context.getEndLocation() - context.getOrigin())) {
//                        int newStart = Math.min((int) context.getOrigin(), interval.getStart());
//                        int newEnd = Math.max((int) context.getEndLocation(), interval.getEnd());
//                        loadedInterval.trimTo(interval.getChr(), newStart, newEnd, interval.getZoom());
//                    }
//
//                    return;
//                }
//            }
//            currentValue.add(interval);
//            while (currentValue.size() > CACHE_SIZE) {
//                currentValue.remove(0);
//            }
//        } else {
//            List<AlignmentInterval> valueList = new LinkedList<AlignmentInterval>();
//            valueList.add(interval);
//            loadedIntervalMap.put(key, valueList);
//        }
    }

    /**
     * TODO -- hacked to get by for now,
     *
     * @return the alignmentRows
     */
    public Map<String, List<AlignmentInterval.Row>> getGroupedAlignments(ReferenceFrame referenceFrame) {
        Collection<AlignmentInterval> loadedIntervals = loadedIntervalMap.get(referenceFrame.getChrName());
        if (loadedIntervals == null) return null;

        Map<String, List<AlignmentInterval.Row>> groupAlignments = new HashMap<String, List<AlignmentInterval.Row>>(loadedIntervals.size());
        for (AlignmentInterval loadedInterval : loadedIntervals) {
            if (loadedInterval.getGroupedAlignments() != null)
                groupAlignments.putAll(loadedInterval.getGroupedAlignments());
        }
        return groupAlignments;
    }

    public int getNLevels() {
        int nLevels = 0;
        for (List<AlignmentInterval> loadedIntervals : loadedIntervalMap.values()) {
            for (AlignmentInterval loadedInterval : loadedIntervals) {
                int intervalNLevels = 0;
                Collection<List<AlignmentInterval.Row>> tmp = loadedInterval.getGroupedAlignments().values();
                for (List<AlignmentInterval.Row> rows : tmp) {
                    intervalNLevels += rows.size();
                }
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
        for (List<AlignmentInterval> loadedIntervals : loadedIntervalMap.values()) {
            for (AlignmentInterval loadedInterval : loadedIntervals) {
                groupCount = Math.max(groupCount, loadedInterval.getGroupCount());
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

    public Collection<AlignmentInterval> getLoadedIntervals() {
        return loadedIntervalMap.getLoadedIntervals();
    }


    public void updatePEStats(AlignmentTrack.RenderOptions renderOptions) {
        if (this.peStats != null) {
            for (PEStats stats : peStats.values()) {
                stats.compute(renderOptions.getMinInsertSizePercentile(), renderOptions.getMaxInsertSizePercentile());
            }
        }
    }

    public boolean isShowSpliceJunctions() {
        return showSpliceJunctions;
    }

    public static class DownsampleOptions {
        private boolean downsample;
        private int sampleWindowSize;
        private int maxReadCount;

        public DownsampleOptions() {
            PreferenceManager prefs = PreferenceManager.getInstance();
            downsample = prefs.getAsBoolean(PreferenceManager.SAM_DOWNSAMPLE_READS);
            sampleWindowSize = prefs.getAsInt(PreferenceManager.SAM_SAMPLING_WINDOW);
            maxReadCount = prefs.getAsInt(PreferenceManager.SAM_SAMPLING_COUNT);
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

