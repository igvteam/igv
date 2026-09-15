package org.igv.alignment.fiberseq;

import org.igv.alignment.AlignmentBlock;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Nucleosomes and methylation sensitive patches (MSPs) called by fibertools, lifted to reference coordinates.
 * <p>
 * fibertools stores them as parallel 0-based read-coordinate arrays: ns/nl (nucleosome starts / lengths) and
 * as/al/aq (MSP starts / lengths / FIRE qualities).  Coordinates are in the original molecule's orientation, so for
 * reverse-strand alignments an interval [s, s+l) corresponds to [L-(s+l), L-s) of the stored read sequence, where L
 * is its length.
 */
public class FiberseqAnnotations {

    /**
     * A 0-based, half-open reference interval.  quality is the FIRE quality (0-255) for MSPs, 0 for nucleosomes.
     */
    public record Interval(int start, int end, int quality) {
        public boolean contains(int position) {
            return position >= start && position < end;
        }
    }

    private final List<Interval> nucleosomes;
    private final List<Interval> msps;

    private FiberseqAnnotations(List<Interval> nucleosomes, List<Interval> msps) {
        this.nucleosomes = nucleosomes;
        this.msps = msps;
    }

    /**
     * Build annotations from raw tag values, or return null if the alignment has neither nucleosome nor MSP tags.
     *
     * @param readLength length of the stored read sequence, including soft clips
     * @param blocks     alignment blocks, whose bases offsets index the stored read sequence
     */
    public static FiberseqAnnotations create(Object ns, Object nl, Object as, Object al, Object aq,
                                             int readLength, boolean isNegativeStrand, AlignmentBlock[] blocks) {
        if (readLength <= 0 || blocks == null) {
            return null;
        }
        List<Interval> nucleosomes = liftIntervals(toIntArray(ns), toIntArray(nl), null, readLength, isNegativeStrand, blocks);
        List<Interval> msps = liftIntervals(toIntArray(as), toIntArray(al), toIntArray(aq), readLength, isNegativeStrand, blocks);
        if (nucleosomes.isEmpty() && msps.isEmpty()) {
            return null;
        }
        return new FiberseqAnnotations(nucleosomes, msps);
    }

    public List<Interval> getNucleosomes() {
        return nucleosomes;
    }

    public List<Interval> getMsps() {
        return msps;
    }

    private static List<Interval> liftIntervals(int[] starts, int[] lengths, int[] qualities, int readLength,
                                                boolean isNegativeStrand, AlignmentBlock[] blocks) {
        if (starts == null || lengths == null || starts.length != lengths.length) {
            return Collections.emptyList();
        }
        List<Interval> intervals = new ArrayList<>(starts.length);
        for (int i = 0; i < starts.length; i++) {
            int start = starts[i];
            int end = start + lengths[i];
            if (isNegativeStrand) {
                int flippedStart = readLength - end;
                end = readLength - start;
                start = flippedStart;
            }
            if (start < 0 || end > readLength || start >= end) {
                continue;
            }
            int[] ref = toReference(blocks, start, end);
            if (ref != null) {
                int quality = qualities != null && i < qualities.length ? qualities[i] : 0;
                intervals.add(new Interval(ref[0], ref[1], quality));
            }
        }
        return intervals;
    }

    /**
     * Lift the read interval [readStart, readEnd) to reference coordinates using the aligned (non soft-clipped)
     * blocks, spanning any deletions or insertions inside it.  Returns null if no aligned base falls in the interval.
     */
    static int[] toReference(AlignmentBlock[] blocks, int readStart, int readEnd) {
        int refStart = -1;
        int refEnd = -1;
        for (AlignmentBlock block : blocks) {
            if (block.isSoftClip()) {
                continue;
            }
            int blockReadStart = block.getBasesOffset();
            int blockReadEnd = blockReadStart + block.getBasesLength();
            if (blockReadEnd <= readStart || blockReadStart >= readEnd) {
                continue;
            }
            if (refStart < 0) {
                refStart = block.getStart() + Math.max(readStart, blockReadStart) - blockReadStart;
            }
            refEnd = block.getStart() + Math.min(readEnd, blockReadEnd) - blockReadStart;
        }
        return refStart < 0 ? null : new int[]{refStart, refEnd};
    }

    /**
     * Convert a B-array tag value to unsigned ints.  fibertools writes ns/nl/as/al as B:I and aq as B:C.
     */
    static int[] toIntArray(Object tagValue) {
        if (tagValue instanceof int[] values) {
            return values;
        } else if (tagValue instanceof short[] values) {
            int[] result = new int[values.length];
            for (int i = 0; i < values.length; i++) result[i] = Short.toUnsignedInt(values[i]);
            return result;
        } else if (tagValue instanceof byte[] values) {
            int[] result = new int[values.length];
            for (int i = 0; i < values.length; i++) result[i] = Byte.toUnsignedInt(values[i]);
            return result;
        }
        return null;
    }
}
