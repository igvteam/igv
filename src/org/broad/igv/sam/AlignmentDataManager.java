/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
package org.broad.igv.sam;

import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.AlignmentTrack.SortOption;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.sam.reader.SamListReader;
import org.broad.igv.track.MultiFileWrapper;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ArrayHeapObjectSorter;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.io.IOException;
import java.util.*;

public class AlignmentDataManager {

    private static Logger log = Logger.getLogger(AlignmentDataManager.class);

    private static final int DEFAULT_DEPTH = 10;

    /**
     * Map of reference frame -> alignment interval
     */
    //TODO -- this is a  potential memory leak, this map needs cleared when the gene list changes
    private HashMap<String, AlignmentInterval> loadedIntervalMap = new HashMap(50);

    HashMap<String, String> chrMappings = new HashMap();
    private boolean isLoading = false;
    private CachingQueryReader reader;
    private CoverageTrack coverageTrack;
    private int maxLevels;


    private boolean viewAsPairs = false;
    private static final int MAX_ROWS = 1000000;
    Map<String, PEStats> peStats;


    public AlignmentDataManager(ResourceLocator locator) throws IOException {

        PreferenceManager prefs = PreferenceManager.getInstance();
        maxLevels = prefs.getAsInt(PreferenceManager.SAM_MAX_LEVELS);

        if (locator.getPath().endsWith(".sam.list")) {
            MultiFileWrapper mfw = MultiFileWrapper.parse(locator);
            reader = new CachingQueryReader(new SamListReader(mfw.getLocators()));
        } else {
            reader = new CachingQueryReader(AlignmentReaderFactory.getReader(locator));
        }
        peStats = new HashMap();
        initChrMap();
    }


    private void initChrMap() {
        Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
        if (genome != null) {
            Set<String> seqNames = reader.getSequenceNames();
            if (seqNames != null) {
                for (String chr : seqNames) {
                    String alias = genome.getChromosomeAlias(chr);
                    chrMappings.put(alias, chr);
                }
            }
        }
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

    public int getMaxLevels() {
        return maxLevels;
    }

    public void setMaxLevels(int maxLevels) {
        clear();
        reader.clearCache();
        this.maxLevels = maxLevels;
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
    public Set<String> getSequenceNames() {
        return reader.getSequenceNames();
    }


    /**
     * Return the loaded interval for the specified frame.  Note this can be null if the interval isn't loaded
     * yet.
     *
     * @param frame
     * @return
     */
    public AlignmentInterval getLoadedInterval(ReferenceFrame frame) {
        return loadedIntervalMap.get(frame.getName());
    }


    public void sortRows(SortOption option, ReferenceFrame referenceFrame, double location, String tag) {
        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        if (loadedInterval != null) {
            loadedInterval.sortRows(option, location, tag);
        }
    }


    public boolean isViewAsPairs() {
        return viewAsPairs;
    }


    public void setViewAsPairs(boolean option, AlignmentTrack.GroupOption groupByOption, String tag) {
        if (option == this.viewAsPairs) {
            return;
        }

        boolean currentPairState = this.viewAsPairs;
        this.viewAsPairs = option;

        for (ReferenceFrame frame : FrameManager.getFrames()) {
            repackAlignments(frame, currentPairState, groupByOption, tag);
        }
    }


    private void repackAlignments(ReferenceFrame referenceFrame, boolean currentPairState,
                                  AlignmentTrack.GroupOption groupByOption, String tag) {
        if (currentPairState == true) {
            AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
            if (loadedInterval == null) {
                return;
            }


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
                    viewAsPairs,
                    groupByOption,
                    tag,
                    MAX_ROWS);

            loadedInterval.setAlignmentRows(tmp);

        } else {
            repackAlignments(referenceFrame, groupByOption, tag);
        }
    }

    /**
     * Repack currently loaded alignments.
     *
     * @param referenceFrame
     */
    public void repackAlignments(ReferenceFrame referenceFrame, AlignmentTrack.GroupOption groupByOption, String tag) {

        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        if (loadedInterval == null) {
            return;
        }

        Iterator<Alignment> iter = loadedInterval.getAlignmentIterator();

        // When repacking keep all currently loaded alignments (don't limit to levels)
        int max = Integer.MAX_VALUE;
        LinkedHashMap<String, List<AlignmentInterval.Row>> alignmentRows = (new AlignmentPacker()).packAlignments(
                iter,
                loadedInterval.getEnd(),
                viewAsPairs,
                groupByOption,
                tag,
                MAX_ROWS);

        loadedInterval.setAlignmentRows(alignmentRows);
    }

    public synchronized LinkedHashMap<String, List<AlignmentInterval.Row>> getGroups(RenderContext context,
                                                                                     AlignmentTrack.GroupOption groupByOption,
                                                                                     String tag) {

        final String genomeId = context.getGenomeId();
        final String chr = context.getChr();
        final int start = (int) context.getOrigin();
        final int end = (int) context.getEndLocation() + 1;

        AlignmentInterval loadedInterval = loadedIntervalMap.get(context.getReferenceFrame().getName());

        // If we've moved out of the loaded interval start a new load.
        if (loadedInterval == null || !loadedInterval.contains(chr, start, end)) {
            loadAlignments(chr, start, end, groupByOption, tag, context);
        }

        // If there is any overlap in the loaded interval and the requested interval return it.
        if (loadedInterval != null && loadedInterval.overlaps(chr, start, end)) {
            return loadedInterval.getGroupedAlignments();
        } else {
            return null;
        }
    }


    public void clear() {
        reader.clearCache();
        loadedIntervalMap.clear();
    }

    public synchronized void loadAlignments(final String chr, final int start, final int end,
                                            final AlignmentTrack.GroupOption groupByOption,
                                            final String tag,
                                            final RenderContext context) {

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

                // Expand start and end to facilitate panning, but by no more than
                // 1 tile size or 8kb  a tile size, whichever is less
                // DON'T expand mitochondria

                int expandLength = Math.min(8000, reader.getTileSize(chr) / 2);
                int intervalStart = start - expandLength;
                int intervalEnd = end + expandLength;

                CloseableIterator<Alignment> iter = null;
                try {

                    String sequence = chrMappings.containsKey(chr) ? chrMappings.get(chr) : chr;

                    List<AlignmentCounts> counts = new ArrayList();

                    iter = reader.query(sequence, intervalStart, intervalEnd, counts, maxLevels, peStats);

                    final AlignmentPacker alignmentPacker = new AlignmentPacker();

                    LinkedHashMap<String, List<AlignmentInterval.Row>> alignmentRows = alignmentPacker.packAlignments(iter,
                            intervalEnd, viewAsPairs, groupByOption, tag, maxLevels);

                    AlignmentInterval loadedInterval = new AlignmentInterval(chr, intervalStart, intervalEnd,
                            alignmentRows, counts);
                    loadedIntervalMap.put(context.getReferenceFrame().getName(), loadedInterval);


                    if (coverageTrack != null) {
                        coverageTrack.rescale(context.getReferenceFrame());
                    }

                    // TODO --- we need to force a repaint of the coverageTrack, which might not be in the same panel
                    if (context.getPanel() != null) {
                        context.getPanel().invalidate();
                        context.getPanel().repaint();
                    }

                } catch (Exception exception) {
                    log.error("Error loading alignments", exception);
                    JOptionPane.showMessageDialog(IGV.getMainFrame(), "Error reading file: " + exception.getMessage());
                } finally {

                    if (iter != null) {
                        iter.close();
                    }
                    isLoading = false;
                }
            }
        };

        LongRunningTask.submit(runnable);


    }


    private boolean isMitochondria(String chr) {
        return chr.equals("M") || chr.equals("chrM") ||
                chr.equals("MT") || chr.equals("chrMT");
    }

    /**
     * TODO -- hacked to get by for now,
     *
     * @return the alignmentRows
     */
    public Map<String, List<AlignmentInterval.Row>> getGroupedAlignments(ReferenceFrame referenceFrame) {
        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        return loadedInterval == null ? null : loadedInterval.getGroupedAlignments();
    }

    public int getNLevels() {
        int nLevels = 0;
        for (AlignmentInterval loadedInterval : loadedIntervalMap.values()) {
            int intervalNLevels = 0;
            Collection<List<AlignmentInterval.Row>> tmp = loadedInterval.getGroupedAlignments().values();
            for (List<AlignmentInterval.Row> rows : tmp) {
                intervalNLevels += rows.size();
            }
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
        for (AlignmentInterval loadedInterval : loadedIntervalMap.values()) {
            groupCount = Math.max(groupCount, loadedInterval.getGroupCount());
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
        return loadedIntervalMap.values();
    }

    public boolean isLoading() {
        return isLoading;
    }

    public void updatePEStats(AlignmentTrack.RenderOptions renderOptions) {
        if (this.peStats != null) {
            for (PEStats stats : peStats.values()) {
                stats.compute(renderOptions.getMinInsertSizePercentile(), renderOptions.getMaxInsertSizePercentile());
            }
        }
    }

    public void updateSpliceJunctions() {
        final Collection<AlignmentInterval> loadedIntervals = getLoadedIntervals();
        if (loadedIntervals != null) {
            for (AlignmentInterval interval : loadedIntervals) {
                interval.resetSpliceJunctions();

            }
        }
    }

}

