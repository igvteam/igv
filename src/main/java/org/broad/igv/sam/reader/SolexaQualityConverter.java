package org.broad.igv.sam.reader;

/**
 * Optimized method for converting Solexa ASCII qualities into Phred scores.
 * Pre-computes all values in order to eliminate repeated computation.
 */
public class SolexaQualityConverter {

    /**
     * This value is added to a Solexa quality score to make it printable ASCII
     */
    public static final int SOLEXA_ADDEND = 64;


    private static SolexaQualityConverter singleton = null;

    public static synchronized SolexaQualityConverter getSingleton() {
        if (singleton == null) {
            singleton = new SolexaQualityConverter();
        }
        return singleton;
    }

    /**
     * Mapping from ASCII value in Gerald export file to phred score
     */
    private final byte[] phredScore = new byte[256];

    private SolexaQualityConverter() {
        for (int i = 0; i < SOLEXA_ADDEND; ++i) {
            phredScore[i] = 0;
        }
        for (int i = SOLEXA_ADDEND; i < phredScore.length; ++i) {
            phredScore[i] = convertSolexaQualityCharToPhredBinary(i);
        }
    }


    /**
     * Converts a solexa character quality into a phred numeric quality.
     */
    private byte convertSolexaQualityCharToPhredBinary(final int solexaQuality) {
        return (byte) Math.round(10d * Math.log10(1d + Math.pow(10d, (solexaQuality - SOLEXA_ADDEND) / 10d)));
    }


    /**
     * Casava 1.3 stores phred-scaled qualities, but non-standard because they have 64 added to them
     * rather than the standard 33.
     *
     * @param solexaQuals qualities are converted in place.
     */
    public void convertSolexa_1_3_QualityCharsToPhredBinary(final byte[] solexaQuals) {
        for (int i = 0; i < solexaQuals.length; ++i) {
            solexaQuals[i] -= SOLEXA_ADDEND;
        }

    }
}
