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
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 */
public class AlignmentInterval extends Locus {

    private static Logger log = Logger.getLogger(AlignmentInterval.class);

    Genome genome;
    private int maxCount = 0;
    private List<AlignmentCounts> counts;
    private LinkedHashMap<String, List<Row>> groupedAlignmentRows;
    private List<SpliceJunctionFeature> spliceJunctions;


    public AlignmentInterval(String chr, int start, int end, LinkedHashMap<String, List<Row>> groupedAlignmentRows,
                             List<AlignmentCounts> counts, List<SpliceJunctionFeature> spliceJunctions) {

        super(chr, start, end);
        this.groupedAlignmentRows = groupedAlignmentRows;
        genome = IGV.getInstance().getGenomeManager().getCurrentGenome();

        //reference = genome.getSequence(chr, start, end);
        this.counts = counts;
        for (AlignmentCounts c : counts) {
            maxCount = Math.max(maxCount, c.getMaxCount());
        }

        this.spliceJunctions = spliceJunctions;

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
    public LinkedHashMap<String, List<Row>> getGroupedAlignments() {
        return groupedAlignmentRows;
    }

    public int getGroupCount() {
        return groupedAlignmentRows == null ? 0 : groupedAlignmentRows.size();
    }

    public void setAlignmentRows(LinkedHashMap<String, List<Row>> alignmentRows) {
        this.groupedAlignmentRows = alignmentRows;
    }


    public void sortRows(AlignmentTrack.SortOption option, ReferenceFrame referenceFrame, String tag) {
        double center = referenceFrame.getCenter();
        sortRows(option, center, tag);
    }


    /**
     * Sort rows group by group
     *
     * @param option
     * @param location
     */
    public void sortRows(AlignmentTrack.SortOption option, double location, String tag) {
        if (groupedAlignmentRows == null) {
            return;
        }

        for (List<AlignmentInterval.Row> alignmentRows : groupedAlignmentRows.values()) {
            for (AlignmentInterval.Row row : alignmentRows) {
                row.updateScore(option, location, this, tag);
            }

            Collections.sort(alignmentRows, new Comparator<Row>() {

                public int compare(AlignmentInterval.Row arg0, AlignmentInterval.Row arg1) {
                    if (arg0.getScore() > arg1.getScore()) {
                        return 1;
                    } else if (arg0.getScore() < arg1.getScore()) {
                        return -1;
                    }
                    return 0;
                }
            });
        }
    }


    public byte getReference(int pos) {
        if (genome == null) {
            return 0;
        }
        return genome.getReference(getChr(), pos);
    }

    public List<AlignmentCounts> getCounts() {
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

    public AlignmentCounts getAlignmentCounts(int pos) {
        for (AlignmentCounts c : counts) {
            if (pos >= c.getStart() && pos < c.getEnd()) {
                return c;
            }
        }
        return null;

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

    public int getDelCount(int pos) {
        for (AlignmentCounts c : counts) {
            if (pos >= c.getStart() && pos < c.getEnd()) {
                return c.getDelCount(pos);
            }
        }
        return 0;
    }

    public int getInsCount(int pos) {
        for (AlignmentCounts c : counts) {
            if (pos >= c.getStart() && pos < c.getEnd()) {
                return c.getInsCount(pos);
            }
        }
        return 0;
    }

    public Iterator<Alignment> getAlignmentIterator() {
        return new AlignmentIterator();
    }

    public List<SpliceJunctionFeature> getSpliceJunctions() throws IOException {
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

        public void updateScore(AlignmentTrack.SortOption option, double center, AlignmentInterval interval, String tag) {

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
                    case FIRST_OF_PAIR_STRAND:
                        Strand strand = centerAlignment.getFirstOfPairStrand();
                        int score = 2;
                        if (strand != Strand.NONE) {
                            score = strand == Strand.NEGATIVE ? 1 : -1;
                        }
                        setScore(score);
                        break;
                    case NUCELOTIDE:
                        byte base = centerAlignment.getBase(adjustedCenter);
                        byte ref = interval.getReference(adjustedCenter);
                        if (base == 0) {      // Base not covered (splice junction)
                            setScore(Integer.MAX_VALUE);
                        } else if (base == 'N' || base == 'n') {
                            setScore(Integer.MAX_VALUE - 2);  // Base is "n"
                        } else if (base == ref) {
                            setScore(Integer.MAX_VALUE - 1);  // Base is reference
                        } else {
                            int count = interval.getCount(adjustedCenter, base);
                            byte phred = centerAlignment.getPhred(adjustedCenter);
                            setScore(-(count + (phred / 100.0f)));
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
                    case MATE_CHR:
                        ReadMate mate = centerAlignment.getMate();
                        if (mate == null) {
                            setScore(Integer.MAX_VALUE);
                        } else {
                            if (mate.getChr().equals(centerAlignment.getChr())) {
                                setScore(Integer.MAX_VALUE - 1);
                            } else {
                                setScore(mate.getChr().hashCode());
                            }
                        }
                        break;
                    case TAG:
                        Object tagValue = centerAlignment.getAttribute(tag);
                        score = tagValue == null ? 0 : tagValue.hashCode();
                        setScore(score);
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

            for (List<AlignmentInterval.Row> alignmentRows : groupedAlignmentRows.values()) {
                for (AlignmentInterval.Row r : alignmentRows) {
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