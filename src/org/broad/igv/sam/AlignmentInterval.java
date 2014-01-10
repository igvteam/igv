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

import org.apache.log4j.Logger;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;

import java.util.*;

/**
 * @author jrobinso
 */
public class AlignmentInterval extends Locus {

    private static Logger log = Logger.getLogger(AlignmentInterval.class);

    Genome genome;
    private AlignmentCounts counts;
    private List<Alignment> alignments;
    private SpliceJunctionHelper spliceJunctionHelper;
    private List<DownsampledInterval> downsampledIntervals;
    private AlignmentTrack.RenderOptions renderOptions;

    AlignmentInterval(AlignmentInterval interval) {
        this(interval.getChr(), interval.getStart(), interval.getEnd(),
                interval.getAlignments(), interval.getCounts(),
                new SpliceJunctionHelper(interval.getSpliceJunctionHelper()), interval.getDownsampledIntervals(), interval.renderOptions);
    }

    public AlignmentInterval(String chr, int start, int end,
                             List<Alignment> alignments,
                             AlignmentCounts counts,
                             SpliceJunctionHelper spliceJunctionHelper,
                             List<DownsampledInterval> downsampledIntervals,
                             AlignmentTrack.RenderOptions renderOptions) {

        super(chr, start, end);
        this.alignments = alignments;
        genome = GenomeManager.getInstance().getCurrentGenome();
        this.counts = counts;

        this.spliceJunctionHelper = spliceJunctionHelper;
        this.downsampledIntervals = downsampledIntervals;
        this.renderOptions = renderOptions;
    }

    static Alignment getFeatureContaining(List<Alignment> features, int right) {

        int leftBounds = 0;
        int rightBounds = features.size() - 1;
        int idx = features.size() / 2;
        int lastIdx = -1;

        while (idx != lastIdx) {
            lastIdx = idx;
            Alignment f = features.get(idx);
            if (f.contains(right)) {
                return f;
            }

            if (f.getStart() > right) {
                rightBounds = idx;
                idx = (leftBounds + idx) / 2;
            } else {
                leftBounds = idx;
                idx = (rightBounds + idx) / 2;

            }

        }
        // Check the extremes
        if (features.get(0).contains(right)) {
            return features.get(0);
        }

        if (features.get(rightBounds).contains(right)) {
            return features.get(rightBounds);
        }

        return null;
    }

    public byte getReference(int pos) {
        if (genome == null) {
            return 0;
        }
        return genome.getReference(getChr(), pos);
    }

    public AlignmentCounts getCounts() {
        return counts;
    }

    /**
     * Return the count of the specified nucleotide
     *
     * @param pos genomic position
     * @param b   nucleotide
     * @return
     */
    public int getCount(int pos, byte b) {
        AlignmentCounts c = counts;
        if (pos >= c.getStart() && pos < c.getEnd()) {
            return c.getCount(pos, b);
        }
        return 0;
    }

    public int getMaxCount(int origin, int end) {
        return counts.getMaxCount(origin, end);
    }

    public int getTotalCount(int pos) {
        AlignmentCounts c = counts;
        if (pos >= c.getStart() && pos < c.getEnd()) {
            return c.getTotalCount(pos);
        }
        return 0;
    }

    public int getDelCount(int pos) {
        AlignmentCounts c = counts;
        if (pos >= c.getStart() && pos < c.getEnd()) {
            return c.getDelCount(pos);
        }
        return 0;
    }

    public List<Alignment> getAlignments(){
        return Collections.unmodifiableList(this.alignments);
    }

    public Iterator<Alignment> getAlignmentIterator() {
        return alignments.iterator();
    }

    public List<DownsampledInterval> getDownsampledIntervals() {
        return downsampledIntervals;
    }

    public SpliceJunctionHelper getSpliceJunctionHelper() {
        return this.spliceJunctionHelper;
    }


    /**
     * An alignment iterator that iterates over packed rows.  Used for
     * repacking.   Using the iterator avoids the need to copy alignments
     * from the rows
     */
    static class AlignmentIterator implements Iterator<Alignment> {

        PriorityQueue<Row> rows;
        Alignment nextAlignment;

        AlignmentIterator(Map<String, List<Row>> groupedAlignmentRows) {
            rows = new PriorityQueue(5, new Comparator<Row>() {

                public int compare(Row o1, Row o2) {
                    return o1.getNextStartPos() - o2.getNextStartPos();
                }
            });

            for (List<Row> alignmentRows : groupedAlignmentRows.values()) {
                for (Row r : alignmentRows) {
                    r.resetIdx();
                    rows.add(r);
                }
            }

            advance();
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
            Row nextRow = null;
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
