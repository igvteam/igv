package org.igv.alignment.fiberseq;

import org.igv.alignment.AlignmentBlock;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

/**
 * Nucleosomes and methylation sensitive patches (MSPs) called by fibertools, lifted to reference coordinates.
 * <p>
 * Two tag encodings are supported.  The molecular annotation tags (Ma:Z and Aq:B:C, or the all-uppercase MA and AQ
 * written by fibertools 0.10-0.12) hold "readLength;type[strand][qualities]:start-length,..." sections with 1-based
 * starts, where each quality letter takes one Aq value per annotation.  FIRE qualities are carried by "fire" entries
 * that duplicate the MSPs called as FIREs.  The legacy tags hold 0-based parallel arrays: ns/nl (nucleosome starts and
 * lengths) and as/al/aq (MSP starts, lengths and FIRE qualities).
 * <p>
 * In both encodings coordinates are in the original molecule's orientation, so for reverse-strand alignments an
 * interval [s, s+l) corresponds to [L-(s+l), L-s) of the stored read sequence, where L is its length.
 */
public class MolecularAnnotations {

    /**
     * A 0-based, half-open reference interval.  quality is the FIRE quality (0-255) for MSPs, 0 for nucleosomes.
     */
    public record Interval(int start, int end, int quality) {
        public boolean contains(int position) {
            return position >= start && position < end;
        }
    }

    /**
     * An annotation covering a position, and the interval carrying it.
     */
    public record Annotation(Type type, Interval interval) {
    }

    /**
     * The annotation an alignment carries at a position, in the order the group and sort options rank them.
     * A FIRE is an MSP with a FIRE quality.
     */
    public enum Type {
        FIRE("FIRE"),
        MSP("MSP"),
        NUCLEOSOME("nucleosome");

        public final String label;

        Type(String label) {
            this.label = label;
        }
    }

    private static final String NUC_TYPE = "nuc";
    private static final String MSP_TYPE = "msp";
    private static final String FIRE_TYPE = "fire";

    private final List<Interval> nucleosomes;
    private final List<Interval> msps;

    private MolecularAnnotations(List<Interval> nucleosomes, List<Interval> msps) {
        this.nucleosomes = nucleosomes;
        this.msps = msps;
    }

    /**
     * True if the tags include fibertools annotations in either encoding.
     *
     * @param tags tag name to value, e.g. {@code SAMRecord::getAttribute}
     */
    public static boolean hasTags(Function<String, Object> tags) {
        return tags.apply("MA") instanceof String || tags.apply("Ma") instanceof String ||
                tags.apply("ns") != null || tags.apply("as") != null;
    }

    /**
     * Build annotations from an alignment's tags, or return null if it has none.  Molecular annotation tags take
     * precedence over legacy tags.  When both spellings are present the uppercase family is used, with its quality
     * tag read in the same spelling, as fibertools does.
     *
     * @param tags      tag name to value, e.g. {@code SAMRecord::getAttribute}
     * @param seqLength length of the stored read sequence including soft clips, 0 if the sequence is absent
     * @param blocks    alignment blocks, whose bases offsets index the stored read sequence
     */
    public static MolecularAnnotations fromTags(Function<String, Object> tags, int seqLength, boolean isNegativeStrand,
                                                AlignmentBlock[] blocks) {
        if (tags.apply("MA") instanceof String ma) {
            return createFromMa(ma, tags.apply("AQ"), seqLength, isNegativeStrand, blocks);
        } else if (tags.apply("Ma") instanceof String ma) {
            return createFromMa(ma, tags.apply("Aq"), seqLength, isNegativeStrand, blocks);
        }
        return create(tags.apply("ns"), tags.apply("nl"), tags.apply("as"), tags.apply("al"), tags.apply("aq"),
                seqLength, isNegativeStrand, blocks);
    }

    public List<Interval> getNucleosomes() {
        return nucleosomes;
    }

    public List<Interval> getMsps() {
        return msps;
    }

    /**
     * The annotation covering a reference position, or null if none does.  MSPs take precedence over nucleosomes
     * should the two overlap.
     */
    public Annotation annotationAt(int position) {
        for (Interval interval : msps) {
            if (interval.contains(position)) {
                return new Annotation(interval.quality() > 0 ? Type.FIRE : Type.MSP, interval);
            }
        }
        for (Interval interval : nucleosomes) {
            if (interval.contains(position)) {
                return new Annotation(Type.NUCLEOSOME, interval);
            }
        }
        return null;
    }

    /**
     * Build annotations from legacy ns/nl/as/al/aq tag values.
     */
    static MolecularAnnotations create(Object ns, Object nl, Object as, Object al, Object aq,
                                       int readLength, boolean isNegativeStrand, AlignmentBlock[] blocks) {
        return build(toIntArray(ns), toIntArray(nl), toIntArray(as), toIntArray(al), toIntArray(aq),
                readLength, isNegativeStrand, blocks);
    }

    /**
     * Build annotations from an Ma tag and its Aq qualities.  Returns null for a malformed tag, or for a stale one
     * whose read length disagrees with the stored sequence (the read was rewritten after tagging).  Without a stored
     * sequence the tag's read length is used to flip reverse-strand intervals.
     */
    static MolecularAnnotations createFromMa(String ma, Object aqTag, int seqLength, boolean isNegativeStrand,
                                             AlignmentBlock[] blocks) {
        int[] aq = aqTag instanceof byte[] ? toIntArray(aqTag) : null;
        List<int[]> nucleosomes = new ArrayList<>();       // {start, length}
        List<int[]> msps = new ArrayList<>();              // {start, length, quality}
        Map<Long, Integer> fireQualities = new HashMap<>();
        int readLength;
        try {
            String[] sections = ma.split(";");
            readLength = Integer.parseInt(sections[0]);
            int aqIndex = 0;
            for (int s = 1; s < sections.length; s++) {
                String section = sections[s];
                if (section.isEmpty()) {
                    continue;
                }
                int colon = section.indexOf(':');
                int strand = colon < 0 ? -1 : strandIndex(section.substring(0, colon));
                if (strand <= 0) {
                    return null;
                }
                String type = section.substring(0, strand);
                String qualitySpec = section.substring(strand + 1, colon);
                if (!qualitySpec.matches("[PQ]*")) {
                    return null;
                }
                int nQualities = qualitySpec.length();
                for (String position : section.substring(colon + 1).split(",")) {
                    if (position.isEmpty()) {
                        continue;
                    }
                    int dash = position.indexOf('-');
                    if (dash < 0) {
                        return null;
                    }
                    int start = Integer.parseInt(position.substring(0, dash)) - 1;
                    int length = Integer.parseInt(position.substring(dash + 1));
                    if (start < 0 || length < 0) {
                        return null;
                    }
                    int quality = 0;
                    if (nQualities > 0) {
                        if (aq == null || aqIndex + nQualities > aq.length) {
                            return null;
                        }
                        // fibertools types carry a single quality; use the first
                        quality = aq[aqIndex];
                        aqIndex += nQualities;
                    }
                    switch (type) {
                        case NUC_TYPE -> nucleosomes.add(new int[]{start, length});
                        case MSP_TYPE -> msps.add(new int[]{start, length, quality});
                        case FIRE_TYPE -> fireQualities.put(intervalKey(start, length), quality);
                        default -> {
                        }
                    }
                }
            }
        } catch (NumberFormatException e) {
            return null;
        }
        if (seqLength > 0 && seqLength != readLength) {
            return null;
        }

        // FIRE entries are exact copies of the MSPs called as FIREs; other MSPs keep their own quality
        int[] mspQualities = new int[msps.size()];
        for (int i = 0; i < mspQualities.length; i++) {
            int[] msp = msps.get(i);
            mspQualities[i] = fireQualities.getOrDefault(intervalKey(msp[0], msp[1]), msp[2]);
        }
        return build(column(nucleosomes, 0), column(nucleosomes, 1), column(msps, 0), column(msps, 1), mspQualities,
                readLength, isNegativeStrand, blocks);
    }

    private static MolecularAnnotations build(int[] nucStarts, int[] nucLengths, int[] mspStarts, int[] mspLengths,
                                              int[] mspQualities, int readLength, boolean isNegativeStrand,
                                              AlignmentBlock[] blocks) {
        if (readLength <= 0 || blocks == null) {
            return null;
        }
        List<Interval> nucleosomes = liftIntervals(nucStarts, nucLengths, null, readLength, isNegativeStrand, blocks);
        List<Interval> msps = liftIntervals(mspStarts, mspLengths, mspQualities, readLength, isNegativeStrand, blocks);
        if (nucleosomes.isEmpty() && msps.isEmpty()) {
            return null;
        }
        return new MolecularAnnotations(nucleosomes, msps);
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
     * Convert a B-array tag value to unsigned ints.  fibertools writes ns/nl/as/al as B:I and aq, Aq as B:C.
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

    private static int strandIndex(String header) {
        for (int i = 0; i < header.length(); i++) {
            char c = header.charAt(i);
            if (c == '+' || c == '-' || c == '.') {
                return i;
            }
        }
        return -1;
    }

    private static long intervalKey(int start, int length) {
        return ((long) start << 32) | (length & 0xffffffffL);
    }

    private static int[] column(List<int[]> rows, int index) {
        int[] values = new int[rows.size()];
        for (int i = 0; i < values.length; i++) {
            values[i] = rows.get(i)[index];
        }
        return values;
    }
}
