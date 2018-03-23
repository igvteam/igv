/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
