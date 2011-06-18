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
