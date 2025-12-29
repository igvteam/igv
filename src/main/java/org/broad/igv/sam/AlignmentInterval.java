package org.broad.igv.sam;

import org.broad.igv.logging.*;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AlignmentInterval extends Locus {

    private static Logger log = LogManager.getLogger(AlignmentInterval.class);
    private  ReferenceFrame referenceFrame;
    Genome genome;
    private AlignmentCounts counts;
    private List<Alignment> alignments;
    private SpliceJunctionHelper spliceJunctionHelper;
    private List<DownsampledInterval> downsampledIntervals;
    private PackedAlignments packedAlignments;

    public AlignmentInterval(String chr, int start, int end,
                             List<Alignment> alignments,
                             AlignmentCounts counts,
                             SpliceJunctionHelper spliceJunctionHelper,
                             List<DownsampledInterval> downsampledIntervals,
                             ReferenceFrame frame) {

        super(chr, start, end);
        this.alignments = alignments;
        genome = GenomeManager.getInstance().getCurrentGenome();
        this.counts = counts;
        this.spliceJunctionHelper = spliceJunctionHelper;
        this.downsampledIntervals = downsampledIntervals;
        this.referenceFrame = frame;
    }

    /**
     * Constructor to create empty interval, that is interval with no alignments
     *
     * @param chr
     * @param start
     * @param end
     */
    public AlignmentInterval(String chr, int start, int end) {
        super(chr, start, end);
        alignments = Collections.EMPTY_LIST;
        counts = EmptyAlignmentCounts.getInstance();
        spliceJunctionHelper = new SpliceJunctionHelper();
        downsampledIntervals = Collections.EMPTY_LIST;
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
     * Sort rows group by group
     *
     * @param option
     * @param location
     * @param tag
     * @param invertSort
     */
    public void sortRows(SortOption option, double location, String tag, boolean invertSort) {

        PackedAlignments packedAlignments = getPackedAlignments();
        if (packedAlignments == null) {
            return;
        }

        final int center = (int) location;
        byte referenceBase = this.getReference(center);
        Comparator<Row> rowComparator = option.getComparator(center, referenceBase, tag, invertSort);

        for (List<Row> alignmentRows : packedAlignments.values()) {
            alignmentRows.sort(rowComparator);
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

    public List<Alignment> getAlignments() {
        return alignments == null ? Collections.<Alignment>emptyList() : Collections.unmodifiableList(alignments);
    }

    public Iterator<Alignment> getAlignmentIterator() {
        return alignments == null ? Collections.<Alignment>emptyList().iterator() : alignments.iterator();
    }

    public List<DownsampledInterval> getDownsampledIntervals() {
        return downsampledIntervals;
    }


    public SpliceJunctionHelper getSpliceJunctionHelper() {
        return this.spliceJunctionHelper;
    }

    public Range getRange() {
        return new Range(getChr(), getStart(), getEnd());
    }

    public void packAlignments(AlignmentTrack.RenderOptions renderOptions, Track.DisplayMode displayMode) {

        final AlignmentPacker alignmentPacker = new AlignmentPacker();
        this.packedAlignments = alignmentPacker.packAlignments(this, renderOptions, this.referenceFrame, displayMode);
    }

    public PackedAlignments getPackedAlignments() {
        return packedAlignments;
    }

    public void dumpAlignments() {
        if (this.alignments != null) this.alignments.clear();
        this.packedAlignments = null;
    }


}
