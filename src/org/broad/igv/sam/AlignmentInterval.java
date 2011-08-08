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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.SequenceManager;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 */
public class AlignmentInterval extends Locus {

    private static Logger log = Logger.getLogger(AlignmentInterval.class);

    private String genomeId;
    private byte[] reference;
    private int maxCount = 0;
    private List<AlignmentCounts> counts;
    private List<AlignmentInterval.Row> alignmentRows;
    private List<SpliceJunctionFeature> spliceJunctions = null;

    public AlignmentInterval(String genomeId, String chr, int start, int end, List<Row> rows, List<AlignmentCounts> counts) {
        super(chr, start, end);
        this.genomeId = genomeId;
        this.alignmentRows = rows;
        reference = SequenceManager.readSequence(this.genomeId, chr, start, end);
        this.counts = counts;
        for (AlignmentCounts c : counts) {
            maxCount = Math.max(maxCount, c.getMaxCount());
        }

        // Force caclulation of splice junctions
        boolean showSpliceJunctionTrack = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_JUNCTION_TRACK);
        if (showSpliceJunctionTrack) {
            try {
                getSpliceJunctions();
            } catch (IOException e) {
                log.error("Error computing splice junctions", e);
            }
        }
    }

    public boolean contains(String genomeId, String chr, int start, int end) {
        return this.genomeId.equals(genomeId) && contains(chr, start, end);
    }


    public boolean overlaps(String genomeId, String chr, int start, int end) {
        return this.genomeId.equals(genomeId) && overlaps(chr, start, end);
    }

    public int getDepth() {
        return alignmentRows.size();
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

    /**
     * The "packed" alignments in this interval
     */
    public List<Row> getAlignmentRows() {
        return alignmentRows;
    }

    public void setAlignmentRows(List<Row> alignmentRows) {
        this.alignmentRows = alignmentRows;
    }


    public void sortRows(AlignmentTrack.SortOption option, ReferenceFrame referenceFrame) {
        double center = referenceFrame.getCenter();
        sortRows(option, center);
    }

    public void sortRows(AlignmentTrack.SortOption option, double location) {
        if (alignmentRows == null) {
            return;
        }

        // TODO -- need the context.  This is here so it will compile       
        for (AlignmentInterval.Row row : alignmentRows) {
            if (option == AlignmentTrack.SortOption.NUCELOTIDE) {

            }
            row.updateScore(option, location, this);
        }

        Collections.sort(alignmentRows, new Comparator<Row>() {

            public int compare(AlignmentInterval.Row arg0, AlignmentInterval.Row arg1) {
                if (arg0.getScore() > arg1.getScore()) {
                    return 1;
                } else if (arg0.getScore() > arg1.getScore()) {
                    return -1;
                }
                return 0;
            }
        });
    }


    public byte getReference(int pos) {
        if (reference == null) {
            return 0;
        }
        int offset = pos - start;
        if (offset < 0 || offset >= reference.length) {
            if (log.isDebugEnabled()) {
                log.debug("Position out of range: " + pos + " (valid range - " + start + "-" + end);
            }
            return 0;
        } else {
            return reference[offset];
        }
    }

    public List<AlignmentCounts> getCounts() {
        return counts;
    }

    /**
     * Return the count of the specified nucleotide
     * @param pos genomic position
     * @param b nucleotide
     * @return
     */
    public int getCount(int pos, byte b) {
        for (AlignmentCounts c : counts) {
            if (pos >= c.getStart() && pos < c.getEnd()) {
                return c.getCount(pos, b);
            }
        }
        return 0;
    }

    public int getMaxCount() {
        return maxCount;
    }

    public int getTotalCount(int pos) {
        for (AlignmentCounts c : counts) {
            if (pos >= c.getStart() && pos < c.getEnd()) {
                return c.getTotalCount(pos);
            }
        }
        return 0;
    }

    public int getNegCount(int pos, byte b) {
        for (AlignmentCounts c : counts) {
            if (pos >= c.getStart() && pos < c.getEnd()) {
                return c.getNegCount(pos, b);
            }
        }
        return 0;
    }

    public int getPosCount(int pos, byte b) {
        for (AlignmentCounts c : counts) {
            if (pos >= c.getStart() && pos < c.getEnd()) {
                return c.getPosCount(pos, b);
            }
        }
        return 0;
    }

    public Iterator<Alignment> getAlignmentIterator() {
        return new AlignmentIterator();
    }

    public List<SpliceJunctionFeature> getSpliceJunctions() throws IOException {
        if (spliceJunctions == null) {
            spliceJunctions = SpliceJunctionHelper.computeFeatures(getAlignmentIterator());
        }
        return spliceJunctions;
    }


    public static class Row {
        int nextIdx;
        private double score = 0;
        List<Alignment> alignments;
        private int start;
        private int lastEnd;

        public Row() {
            nextIdx = 0;
            this.alignments = new ArrayList(100);
        }

        public void addAlignment(Alignment alignment) {
            if (alignments.isEmpty()) {
                this.start = alignment.getStart();
            }
            alignments.add(alignment);
            lastEnd = alignment.getEnd();

        }

        public void updateScore(AlignmentTrack.SortOption option, double center, AlignmentInterval interval) {

            int adjustedCenter = (int) center;
            Alignment centerAlignment = getFeatureContaining(alignments, adjustedCenter);
            if (centerAlignment == null) {
                setScore(Double.MAX_VALUE);
            } else {
                switch (option) {
                    case START:
                        setScore(centerAlignment.getStart());
                        break;
                    case STRAND:
                        setScore(centerAlignment.isNegativeStrand() ? -1 : 1);
                        break;
                    case NUCELOTIDE:
                        byte base = centerAlignment.getBase(adjustedCenter);
                        byte ref = interval.getReference(adjustedCenter);
                        if (base == 'N' || base == 'n') {
                            setScore(Integer.MAX_VALUE - 1);
                        } else if (base == ref) {
                            setScore(Integer.MAX_VALUE);
                        } else {
                            int count = interval.getCount(adjustedCenter, base);
                            byte phred = centerAlignment.getPhred(adjustedCenter);
                            float score = -(count + (phred / 100.0f));
                            setScore(score);
                        }
                        break;
                    case QUALITY:
                        setScore(-centerAlignment.getMappingQuality());
                        break;
                    case SAMPLE:
                        String sample = centerAlignment.getSample();
                        score = sample == null ? 0 : sample.hashCode();
                        setScore(score);
                        break;
                    case READ_GROUP:
                        String readGroup = centerAlignment.getReadGroup();
                        score = readGroup == null ? 0 : readGroup.hashCode();
                        setScore(score);
                        break;
                    case INSERT_SIZE:
                        setScore(-Math.abs(centerAlignment.getInferredInsertSize()));
                        break;
                }
            }
        }


        // Used for iterating over all alignments, e.g. for packing

        public Alignment nextAlignment() {
            if (nextIdx < alignments.size()) {
                Alignment tmp = alignments.get(nextIdx);
                nextIdx++;
                return tmp;
            } else {
                return null;
            }
        }

        public int getNextStartPos() {
            if (nextIdx < alignments.size()) {
                return alignments.get(nextIdx).getStart();
            } else {
                return Integer.MAX_VALUE;
            }
        }

        public boolean hasNext() {
            return nextIdx < alignments.size();
        }

        public void resetIdx() {
            nextIdx = 0;
        }

        /**
         * @return the score
         */
        public double getScore() {
            return score;
        }

        /**
         * @param score the score to set
         */
        public void setScore(double score) {
            this.score = score;
        }

        public int getStart() {
            return start;
        }

        public int getLastEnd() {
            return lastEnd;
        }

    } // end class row


    /**
     * An alignment iterator that iterates over packed rows.  Used for
     * "repacking".   Using the iterator avoids the need to copy alignments
     * from the rows
     */
    class AlignmentIterator implements Iterator<Alignment> {

        PriorityQueue<AlignmentInterval.Row> rows;
        Alignment nextAlignment;

        AlignmentIterator() {
            rows = new PriorityQueue(5, new Comparator<AlignmentInterval.Row>() {

                public int compare(AlignmentInterval.Row o1, AlignmentInterval.Row o2) {
                    return o1.getNextStartPos() - o2.getNextStartPos();
                }
            });

            for (AlignmentInterval.Row r : alignmentRows) {
                r.resetIdx();
                rows.add(r);
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