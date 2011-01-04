/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.AlignmentTrack.SortOption;
import org.broad.igv.sam.reader.CachingQueryReader;
import org.broad.igv.sam.reader.SamListReader;
import org.broad.igv.sam.reader.SamQueryReaderFactory;
import org.broad.igv.track.MultiFileWrapper;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.io.IOException;
import java.util.*;

public class AlignmentDataManager {

    private static Logger log = Logger.getLogger(AlignmentDataManager.class);

    private static final int DEFAULT_DEPTH = 10;

    //TODO -- this is a  memory leak, this map needs cleared when the gene list changes
    private HashMap<String, AlignmentInterval> loadedIntervalMap = new HashMap(50);
    private HashMap<String, String> chrMappings = new HashMap();
    private boolean isLoading = false;
    private CachingQueryReader reader;
    private CoverageTrack coverageTrack;
    private int maxLevels;


    public AlignmentDataManager(ResourceLocator locator) throws IOException {

        PreferenceManager prefs = PreferenceManager.getInstance();
        maxLevels = prefs.getAsInt(PreferenceManager.SAM_MAX_LEVELS);

        if (locator.getPath().endsWith(".sam.list")) {
            MultiFileWrapper mfw = MultiFileWrapper.parse(locator);
            reader = new CachingQueryReader(new SamListReader(mfw.getLocators()));
        } else {
            reader = new CachingQueryReader(SamQueryReaderFactory.getReader(locator));
        }
        initChrMap();
    }


    private void initChrMap() {
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
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


    public int getMaxDepth(ReferenceFrame referenceFrame) {
        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        return loadedInterval == null ? DEFAULT_DEPTH :
                (loadedInterval.getDepth() == 0 ? DEFAULT_DEPTH : loadedInterval.getDepth());
    }

    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
    }

    /**
     * Return the loaded interval for the specified frame.  Note this can be null if the interval isn't loaded
     * yet.
     *
     * @param frame
     * @return
     */
    public AlignmentInterval getLoadedInterval(ReferenceFrame frame) {
        return loadedIntervalMap.get(frame);
    }

    /**
     * Return the loaded interval for the RenderContext.  This method forces a load if the interval isn't present.
     * It is provided to work aroud problems with batch mode image generation.
     *
     * @param context
     * @return
     */
    public AlignmentInterval getLoadedInterval(RenderContext context) {
        ReferenceFrame frame = context.getReferenceFrame();
        if (!loadedIntervalMap.containsKey(frame.getName())) {
            // IF in batch mode force a load of the interval if its missing
            if (Globals.batch) {
                int start = Math.max(0, (int) context.getOrigin() - 100);
                int end = (int) context.getEndLocation() + 100;
                loadAlignments(context.getGenomeId(), context.getChr(), start, end, context);
            }
        }

        return loadedIntervalMap.get(frame.getName());

    }


    /**
     * Sort alignment rows such that alignments that intersect from the
     * center appear left to right by start position
     */
    public void sortRows(SortOption option, ReferenceFrame referenceFrame) {
        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        if (loadedInterval != null) {
            loadedInterval.sortRows(option, referenceFrame);
        }
    }

    public void sortRows(SortOption option, ReferenceFrame referenceFrame, double location) {
        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        if (loadedInterval != null) {
            loadedInterval.sortRows(option, location);
        }
    }


    /**
     * Repack currently loaded alignments.
     *
     * @param referenceFrame
     */
    public void repackAlignments(ReferenceFrame referenceFrame) {
        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        if (loadedInterval == null) {
            return;
        }
        RowIterator iter = new RowIterator(referenceFrame);

        List<AlignmentInterval.Row> alignmentRows = AlignmentLoader.loadAndPackAlignments(
                iter,
                maxLevels,
                loadedInterval.getEnd());

        loadedInterval.setAlignmentRows(alignmentRows);
    }


    public synchronized List<AlignmentInterval.Row> getAlignmentRows(RenderContext context) {

        final String genomeId = context.getGenomeId();
        final String chr = context.getChr();
        final int start = (int) context.getOrigin();
        final int end = (int) context.getEndLocation() + 1;

        AlignmentInterval loadedInterval = loadedIntervalMap.get(context.getReferenceFrame().getName());

        // If we've moved out of the loaded interval start a new load.
        if (loadedInterval == null || !loadedInterval.contains(genomeId, chr, start, end)) {
            log.debug("Loading alignments: " + chr + ":" + start + "-" + end);
            loadAlignments(genomeId, chr, start, end, context);
        }

        // If there is any overlap in the loaded interval and the requested interval return it.
        if (loadedInterval != null && loadedInterval.overlaps(genomeId, chr, start, end)) {
            return loadedInterval.getAlignmentRows();
        } else {
            return null;
        }
    }

    public void clear() {
        reader.clearCache();
        loadedIntervalMap.clear();
    }

    public void loadAlignments(final String genomeId, final String chr, final int start, final int end, final RenderContext context) {

        if (isLoading || chr.equals(Globals.CHR_ALL)) {
            return;
        }

        log.debug("Load alignments.  isLoading=" + isLoading);
        isLoading = true;
        NamedRunnable runnable = new NamedRunnable() {

            public String getName() {
                return "preloadData";
            }

            public void run() {

                 final PreferenceManager prefs = PreferenceManager.getInstance();
                final int qualityThreshold = prefs.getAsInt(PreferenceManager.SAM_QUALITY_THRESHOLD);


                // Expand start and end to facilitate panning, but by no more than
                // 1 screen or 8kb, whichever is less
                // DON'T expand mitochondria

                int expandLength = reader.getTileSize(chr) / 2;
                int intervalStart = Math.max(0, start - expandLength);
                int intervalEnd = end + expandLength;
                CloseableIterator<Alignment> iter = null;
                try {

                    String sequence = chrMappings.containsKey(chr) ? chrMappings.get(chr) : chr;

                    List<CachingQueryReader.AlignmentCounts> counts = new ArrayList();

                    iter = reader.query(sequence, intervalStart, intervalEnd, counts);

                    List<AlignmentInterval.Row> alignmentRows = AlignmentLoader.loadAndPackAlignments(iter, maxLevels,
                            intervalEnd);

                    AlignmentInterval loadedInterval = new AlignmentInterval(genomeId, chr, intervalStart, intervalEnd,
                            alignmentRows, counts);
                    loadedIntervalMap.put(context.getReferenceFrame().getName(), loadedInterval);


                    if (coverageTrack != null) {
                        coverageTrack.rescale(context.getReferenceFrame());
                    }

                    // TODO --- we need to force a repaint of the coverageTrack, which might not be in the same panel
                    if (context.getPanel() != null) context.getPanel().repaint();

                    //TODO -- this has to be done after every load in every panel.  Centralize this somewhere?  Have
                    //TODO --  a "DataLoadRunnable"?
                    IGVMainFrame.getInstance().layoutMainPanel();


                } catch (Exception exception) {
                    log.error("Error loading alignments", exception);
                    JOptionPane.showMessageDialog(IGVMainFrame.getInstance(), "Error reading file: " + exception.getMessage());
                } finally {
                    isLoading = false;
                    if (iter != null) {
                        iter.close();
                    }
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
    public List<AlignmentInterval.Row> getAlignmentRows(ReferenceFrame referenceFrame) {
        AlignmentInterval loadedInterval = loadedIntervalMap.get(referenceFrame.getName());
        return loadedInterval == null ? null : loadedInterval.getAlignmentRows();
    }

    public int getNLevels() {
        int nLevels = 0;
        for (AlignmentInterval loadedInterval : loadedIntervalMap.values()) {
            nLevels = Math.max(nLevels, loadedInterval.getAlignmentRows().size());
        }
        return nLevels;
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


    /**
     * An alignment iterator that iterates over packed rows.  Used for
     * "repacking".   Using the iterator avoids the need to copy alignments
     * from the rows
     */
    class RowIterator implements CloseableIterator<Alignment> {

        PriorityQueue<AlignmentInterval.Row> rows;
        Alignment nextAlignment;

        RowIterator(ReferenceFrame referenceFrame) {
            rows = new PriorityQueue(5, new Comparator<AlignmentInterval.Row>() {

                public int compare(AlignmentInterval.Row o1, AlignmentInterval.Row o2) {
                    return o1.getNextStartPos() - o2.getNextStartPos();
                }
            });

            for (AlignmentInterval.Row r : getAlignmentRows(referenceFrame)) {
                r.resetIdx();
                rows.add(r);
            }

            advance();
        }

        public void close() {
            // Ignored
        }

        public boolean hasNext() {
            return nextAlignment != null;
        }

        public Alignment next() {
            Alignment tmp = nextAlignment;
            if (tmp != null) {
                advance();
            }
            return tmp;
        }

        private void advance() {

            nextAlignment = null;
            AlignmentInterval.Row nextRow = null;
            while (nextAlignment == null && !rows.isEmpty()) {
                while ((nextRow = rows.poll()) != null) {
                    if (nextRow.hasNext()) {
                        nextAlignment = nextRow.nextAlignment();
                        break;
                    }
                }
            }
            if (nextRow != null && nextAlignment != null) {
                rows.add(nextRow);
            }
        }

        public void remove() {
            // ignore
        }
    }

}

