package org.igv.sam;


import htsjdk.samtools.util.CloseableIterator;
import org.igv.Globals;
import org.igv.sam.reader.AlignmentReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 12/6/11
 */
public class AlignmentUtils {

    public static final byte a = 'a';
    public static final byte c = 'c';
    public static final byte g = 'g';
    public static final byte t = 't';
    public static final byte A = 'A';
    public static final byte C = 'C';
    public static final byte G = 'G';
    public static final byte T = 'T';

    /**
     * Return true if the two bases can be considered a match.  The comparison is case-insensitive, and
     * considers ambiguity codes in the reference.
     *
     * @param refbase
     * @param readbase
     * @return
     */
    public static boolean compareBases(byte refbase, byte readbase) {

        if (readbase == 61) {
            return true;  // By definition, 61 is "equals"
        }
        // Force both bases to upper case
        if (refbase > 90) {
            refbase = (byte) (refbase - 32);
        }
        if (readbase > 90) {
            readbase = (byte) (readbase - 32);
        }
        if (refbase == readbase) {
            return true;
        }
        switch (refbase) {
            case 'N':
                return true; // Everything matches 'N'
            case 'U':
                return readbase == 'T';
            case 'M':
                return readbase == 'A' || readbase == 'C';
            case 'R':
                return readbase == 'A' || readbase == 'G';
            case 'W':
                return readbase == 'A' || readbase == 'T';
            case 'S':
                return readbase == 'C' || readbase == 'G';
            case 'Y':
                return readbase == 'C' || readbase == 'T';
            case 'K':
                return readbase == 'G' || readbase == 'T';
            case 'V':
                return readbase == 'A' || readbase == 'C' || readbase == 'G';
            case 'H':
                return readbase == 'A' || readbase == 'C' || readbase == 'T';
            case 'D':
                return readbase == 'A' || readbase == 'G' || readbase == 'T';
            case 'B':
                return readbase == 'C' || readbase == 'G' || readbase == 'T';
            default:
                return false;
        }
    }

    /**
     * Check whether there is a mismatch between {@code reference[idx]} and {@code read[idx]},
     * guarding against {@code reference} being too short.
     * If we do not have a valid reference we assume a match, that is, NOT a misMatch.
     * <p>
     * Note '=' means indicates a match by definition
     *
     * @param reference
     * @param read
     * @param isSoftClipped
     * @param idx
     * @return
     */
    static boolean isMisMatch(byte[] reference, ByteSubarray read, boolean isSoftClipped, int idx) {
        if (reference == null) return false;
        boolean misMatch = false;
        if (isSoftClipped) {
            // Goby will return '=' characters when the soft-clip happens to match the reference.
            // It could actually be useful to see which part of the soft clipped bases match, to help detect
            // cases when an aligner clipped too much.
            final byte readbase = read.getByte(idx);
            misMatch = readbase != '=';  // mismatch, except when the soft-clip has an '=' base.
        } else {
            final int referenceLength = reference.length;
            final byte refbase = idx < referenceLength ? reference[idx] : 0;
            final byte readbase = read.getByte(idx);
            misMatch = readbase != '=' &&
                    idx < referenceLength &&
                    refbase != 0 &&
                    !AlignmentUtils.compareBases(refbase, readbase);
        }
        return misMatch;
    }

    /**
     * Reverses and complements a copy of the original array
     */
    public static byte[] reverseComplementCopy(final byte[] bases) {
        final int lastIndex = bases.length - 1;
        byte[] out = new byte[bases.length];
        int i;
        for (i = 0; i <= lastIndex; i++) {
            out[lastIndex - i] = complement(bases[i]);
        }
        return out;
    }

    /**
     * Reverses and complements the bases in place.
     */
    public static void reverseComplement(final byte[] bases) {
        final int lastIndex = bases.length - 1;

        int i, j;
        for (i = 0, j = lastIndex; i < j; ++i, --j) {
            final byte tmp = complement(bases[i]);
            bases[i] = complement(bases[j]);
            bases[j] = tmp;
        }
        if (bases.length % 2 == 1) {
            bases[i] = complement(bases[i]);
        }
    }

    /**
     * Reverses and complements a copy of the original array
     */
    public static ByteSubarray reverseComplementCopy(ByteSubarray bases) {
        final int lastIndex = bases.length - 1;
        byte[] out = new byte[bases.length];
        int i;
        for (i = 0; i <= lastIndex; i++) {
            out[lastIndex - i] = complement(bases.getByte(i));
        }
        return new ByteSubarray(out, 0, out.length, bases.fillByte);
    }

    /**
     * Returns the complement of a single byte.
     */
    public static final byte complement(final byte b) {
        switch (b) {
            case a:
                return t;
            case c:
                return g;
            case g:
                return c;
            case t:
                return a;
            case A:
                return T;
            case C:
                return G;
            case G:
                return C;
            case T:
                return A;
            default:
                return b;
        }
    }

    /**
     * Calculate the reverse complement of the specified sequence
     * (Stolen from Reseq)
     *
     * @param sequenceData
     * @return reverse complement
     */
    public static String reverseComplement(final String sequenceData) {
        final byte[] bases = htsjdk.samtools.util.StringUtil.stringToBytes(sequenceData);
        reverseComplement(bases);
        return htsjdk.samtools.util.StringUtil.bytesToString(bases);
    }


    public static AlignmentBlockImpl buildAlignmentBlock(char operator,
                                                         byte[] readBases,
                                                         byte[] readBaseQualities,
                                                         int blockStart,
                                                         int fromIdx,
                                                         int nBases) {

        AlignmentBlockImpl block = new AlignmentBlockImpl(blockStart, readBases, readBaseQualities, fromIdx, nBases, operator);

        return block;
    }


    public static List<Alignment> firstAlignments(AlignmentReader reader, int count) throws IOException {
        
        List<Alignment> alignments = new ArrayList<>();
        CloseableIterator<Alignment> iter = null;

        try {
            iter = reader.iterator();
            while (iter.hasNext() && alignments.size() < count) {
                alignments.add(iter.next());
            }
            return alignments;
        } finally {
            iter.close();
        }
    }


    static int[] emptyArray = {};

}
