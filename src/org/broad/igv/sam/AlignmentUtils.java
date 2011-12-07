package org.broad.igv.sam;

import org.broad.igv.sam.reader.GeraldParser;

/**
 * @author Jim Robinson
 * @date 12/6/11
 */
public class AlignmentUtils {
    /**
     * Return true if the two bases can be considered a match.  The comparison is case-insentive, and
     * considers ambiguity codes in the reference.
     *
     * @param refbase
     * @param readbase
     * @return
     */
    public static boolean compareBases(byte refbase, byte readbase) {

        if(readbase == 61) {
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
                return refbase == readbase;
        }
    }

    /**
     * Reverses and complements a copy of the original array
     */
    public static byte[] reverseComplementCopy(final byte[] bases) {
    	final int lastIndex = bases.length - 1;
    	byte[] out = new byte[bases.length];
    	int i;
    	for (i=0; i <= lastIndex; i++)
    	{
    		out[lastIndex-i] = complement(bases[i]);
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
     * Returns the complement of a single byte.
     */
    public static final byte complement(final byte b) {
        switch (b) {
            case GeraldParser.a:
                return GeraldParser.t;
            case GeraldParser.c:
                return GeraldParser.g;
            case GeraldParser.g:
                return GeraldParser.c;
            case GeraldParser.t:
                return GeraldParser.a;
            case GeraldParser.A:
                return GeraldParser.T;
            case GeraldParser.C:
                return GeraldParser.G;
            case GeraldParser.G:
                return GeraldParser.C;
            case GeraldParser.T:
                return GeraldParser.A;
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
        final byte[] bases = net.sf.samtools.util.StringUtil.stringToBytes(sequenceData);
        reverseComplement(bases);
        return net.sf.samtools.util.StringUtil.bytesToString(bases);
    }
}
