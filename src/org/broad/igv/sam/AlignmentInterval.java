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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import com.google.common.base.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.data.Interval;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.SpliceJunctionFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.collections.CollUtils;
import org.broad.tribble.Feature;

import java.util.*;

/**
 * @author jrobinso
 */
public class AlignmentInterval extends Locus implements Interval {

    private static Logger log = Logger.getLogger(AlignmentInterval.class);

    Genome genome;
    private int maxCount = 0;
    private AlignmentCounts counts;
    private LinkedHashMap<String, List<Row>> groupedAlignmentRows;  // The alignments
    private List<SpliceJunctionFeature> spliceJunctions;
    private List<DownsampledInterval> downsampledIntervals;
    private AlignmentTrack.RenderOptions renderOptions;

    public AlignmentInterval(String chr, int start, int end,
                             LinkedHashMap<String, List<Row>> groupedAlignmentRows,
                             AlignmentCounts counts,
                             List<SpliceJunctionFeature> spliceJunctions,
                             List<DownsampledInterval> downsampledIntervals,
                             AlignmentTrack.RenderOptions renderOptions) {

        super(chr, start, end);
        this.groupedAlignmentRows = groupedAlignmentRows;
        genome = GenomeManager.getInstance().getCurrentGenome();

        //reference = genome.getSequence(chr, start, end);
        this.counts = counts;
        this.maxCount = counts.getMaxCount();

        this.spliceJunctions = spliceJunctions;
        this.downsampledIntervals = downsampledIntervals;
        this.renderOptions = renderOptions;
    }

    static AlignmentInterval emptyAlignmentInterval(String chr, int start, int end) {
        return new AlignmentInterval(chr, start, end, null, null, null, null, null);
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

    public void setAlignmentRows(LinkedHashMap<String, List<Row>> alignmentRows, AlignmentTrack.RenderOptions renderOptions) {
        this.groupedAlignmentRows = alignmentRows;
        this.renderOptions = renderOptions;
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

            Collections.sort(alignmentRows);
        }
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

    public int getMaxCount() {
        return maxCount;
    }

    public AlignmentCounts getAlignmentCounts(int pos) {
        AlignmentCounts c = counts;
        if (pos >= c.getStart() && pos < c.getEnd()) {
            return c;
        }
        return null;

    }

    public int getTotalCount(int pos) {
        AlignmentCounts c = counts;
        if (pos >= c.getStart() && pos < c.getEnd()) {
            return c.getTotalCount(pos);
        }
        return 0;
    }

    public int getNegCount(int pos, byte b) {
        AlignmentCounts c = counts;
        if (pos >= c.getStart() && pos < c.getEnd()) {
            return c.getNegCount(pos, b);
        }
        return 0;
    }

    public int getPosCount(int pos, byte b) {
        AlignmentCounts c = counts;
        if (pos >= c.getStart() && pos < c.getEnd()) {
            return c.getPosCount(pos, b);
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

    public int getInsCount(int pos) {
        AlignmentCounts c = counts;
        if (pos >= c.getStart() && pos < c.getEnd()) {
            return c.getInsCount(pos);
        }
        return 0;
    }

    public Iterator<Alignment> getAlignmentIterator() {
        return new AlignmentIterator();
    }

    public List<SpliceJunctionFeature> getSpliceJunctions() {
        return spliceJunctions;
    }

    public List<DownsampledInterval> getDownsampledIntervals() {
        return downsampledIntervals;
    }

    @Override
    public boolean contains(String chr, int start, int end, int zoom) {
        return super.contains(chr, start, end);
    }

    @Override
    public boolean overlaps(String chr, int start, int end, int zoom) {
        return super.overlaps(chr, start, end);
    }

    @Override
    public boolean canMerge(Interval i) {
        if (!super.overlaps(i.getChr(), i.getStart(), i.getEnd())
                || !(i instanceof AlignmentInterval)) {
            return false;
        }
        return getCounts().getClass().equals(((AlignmentInterval) i).getCounts().getClass());
    }

    @Override
    public boolean merge(Interval i) {
        boolean canMerge = this.canMerge(i);
        if (!canMerge) return false;

        AlignmentInterval other = (AlignmentInterval) i;

        List<Alignment> allAlignments = (List<Alignment>) FeatureUtils.combineSortedFeatureListsNoDups(
                getAlignmentIterator(), other.getAlignmentIterator(), start, end);

        this.counts = counts.merge(other.getCounts(), renderOptions.bisulfiteContext);
        this.spliceJunctions = FeatureUtils.combineSortedFeatureListsNoDups(this.spliceJunctions, other.getSpliceJunctions(), start, end);
        this.downsampledIntervals = FeatureUtils.combineSortedFeatureListsNoDups(this.downsampledIntervals, other.getDownsampledIntervals(), start, end);

        //This must be done AFTER calling combineSortedFeatureListsNoDups for the last time,
        //because we rely on the original start/end
        this.start = Math.min(getStart(), i.getStart());
        this.end = Math.max(getEnd(), i.getEnd());
        this.maxCount = Math.max(this.getMaxCount(), other.getMaxCount());


        AlignmentPacker packer = new AlignmentPacker();
        this.groupedAlignmentRows = packer.packAlignments(allAlignments.iterator(), this.end, renderOptions);


        return true;
    }

    /**
     * Remove data outside the specified interval.
     * This interval must contain the specified interval.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @return true if any trimming was performed, false if not
     */
    public boolean trimTo(final String chr, final int start, final int end, int zoom) {
        boolean toRet = false;
        if (!this.contains(chr, start, end, zoom)) {
            return false;
        }

        for (String groupKey : groupedAlignmentRows.keySet()) {
            for (AlignmentInterval.Row row : groupedAlignmentRows.get(groupKey)) {
                Iterator<Alignment> alIter = row.alignments.iterator();
                Alignment al;
                while (alIter.hasNext()) {
                    al = alIter.next();
                    if (al.getEnd() < start || al.getStart() > end) {
                        alIter.remove();
                        toRet |= true;
                    }
                }

            }
        }

        Predicate<Feature> overlapPredicate = FeatureUtils.getOverlapPredicate(chr, start, end);

        if (getDownsampledIntervals() != null) {
            downsampledIntervals = CollUtils.filter(this.getDownsampledIntervals(), overlapPredicate);
        }
        if (getSpliceJunctions() != null) {
            spliceJunctions = CollUtils.filter(this.getSpliceJunctions(), overlapPredicate);
        }

        this.start = Math.max(this.start, start);
        this.end = Math.min(this.end, end);


        return toRet;
    }

    private List addToListNoDups(List self, List other) {
        if (self == null) self = new ArrayList();
        if (other != null) {
            Set selfSet = new HashSet(self);
            selfSet.addAll(other);

            self = new ArrayList(selfSet);
            FeatureUtils.sortFeatureList(self);
        }
        return self;
    }

    /**
     * AlignmentInterval data is independent of zoom
     *
     * @return
     */
    @Override
    public int getZoom() {
        return -1;
    }

    public static class Row implements Comparable<Row> {
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

        public void updateScore(AlignmentTrack.SortOption option, Locus locus, AlignmentInterval interval, String tag) {
            double mean = 0;
            //double sd = 0;
            int number = 0;
            for (int center = locus.getStart(); center < locus.getEnd(); center++) {
                double value = calculateScore(option, center, interval, tag);

                mean = number * (mean / (number + 1)) + (value / (number + 1));

                number++;
            }
            setScore(mean);
        }

        public void updateScore(AlignmentTrack.SortOption option, double center, AlignmentInterval interval, String tag) {
            setScore(calculateScore(option, center, interval, tag));
        }


        public double calculateScore(AlignmentTrack.SortOption option, double center, AlignmentInterval interval, String tag) {

            int adjustedCenter = (int) center;
            Alignment centerAlignment = getFeatureContaining(alignments, adjustedCenter);
            if (centerAlignment == null) {
                return Integer.MAX_VALUE;
            } else {
                switch (option) {
                    case START:
                        return centerAlignment.getStart();
                    case STRAND:
                        return centerAlignment.isNegativeStrand() ? -1 : 1;
                    case FIRST_OF_PAIR_STRAND:
                        Strand strand = centerAlignment.getFirstOfPairStrand();
                        int score = 2;
                        if (strand != Strand.NONE) {
                            score = strand == Strand.NEGATIVE ? 1 : -1;
                        }
                        return score;
                    case NUCELOTIDE:
                        byte base = centerAlignment.getBase(adjustedCenter);
                        byte ref = interval.getReference(adjustedCenter);

                        if (base == 'N' || base == 'n') {
                            return 2;  // Base is "n"
                        } else if (base == ref) {
                            return 3;  // Base is reference
                        } else {
                            //If base is 0, base not covered (splice junction) or is deletion
                            if (base == 0) {
                                int delCount = interval.getDelCount(adjustedCenter);
                                if (delCount > 0) {
                                    return -delCount;
                                } else {
                                    //Base not covered, NOT a deletion
                                    return 1;
                                }
                            } else {
                                int count = interval.getCount(adjustedCenter, base);
                                byte phred = centerAlignment.getPhred(adjustedCenter);
                                return -(count + (phred / 100.0f));
                            }
                        }
                    case QUALITY:
                        return -centerAlignment.getMappingQuality();
                    case SAMPLE:
                        String sample = centerAlignment.getSample();
                        score = sample == null ? 0 : sample.hashCode();
                        return score;
                    case READ_GROUP:
                        String readGroup = centerAlignment.getReadGroup();
                        score = readGroup == null ? 0 : readGroup.hashCode();
                        return score;
                    case INSERT_SIZE:
                        return -Math.abs(centerAlignment.getInferredInsertSize());
                    case MATE_CHR:
                        ReadMate mate = centerAlignment.getMate();
                        if (mate == null) {
                            return Integer.MAX_VALUE;
                        } else {
                            if (mate.getChr().equals(centerAlignment.getChr())) {
                                return Integer.MAX_VALUE - 1;
                            } else {
                                return mate.getChr().hashCode();
                            }
                        }
                    case TAG:
                        Object tagValue = centerAlignment.getAttribute(tag);
                        score = tagValue == null ? 0 : tagValue.hashCode();
                        return score;
                    default:
                        return Integer.MAX_VALUE;
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

        @Override
        public int compareTo(Row o) {
            return (int) Math.signum(getScore() - o.getScore());
        }

//        @Override
//        public boolean equals(Object object){
//            if(!(object instanceof Row)){
//                return false;
//            }
//            Row other = (Row) object;
//            boolean equals = this.getStart() == other.getStart();
//            equals &= this.getLastEnd() == other.getLastEnd();
//            equals &= this.getScore() == other.getScore();
//
//            return equals;
//
//        }
//
//        @Override
//        public int hashCode(){
//            int score = (int) getScore();
//            score = score != 0 ? score : 1;
//            return (getStart() * getLastEnd() * score);
//        }

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
