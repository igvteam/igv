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
package org.broad.igv.maf;

//~--- non-JDK imports --------------------------------------------------------

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class MAFTile {

    int start;
    int end;
    MASequence refSeq;
    int[] gapAdjustedIdx;
    List<Gap> gaps = new ArrayList();
    /**
     * Map of species id -> aligned sequence
     */
    Map<String, MASequence> alignedSequences;

    public MAFTile() {
    }

    // Constructor for an empty maf tile (no alignments)
    public MAFTile(int start, int end) {
        this.start = start;
        this.end = end;
        this.alignedSequences = null;
        this.gapAdjustedIdx = null;

    }

    public MAFTile(int start, int end, Map<String, String> bases, int[] gapAdjustedCoordinates) {
        this(start, end, bases, gapAdjustedCoordinates, "hg18");
    }

    public MAFTile(int start, int end, Map<String, String> bases, int[] gapAdjustedCoordinates, String refSeqId) {
        this.start = start;
        this.end = end;

        this.gapAdjustedIdx = gapAdjustedCoordinates;

        int expectedGapAdjustedIdx = 0;
        for (int i = start; i < end; i++) {
            int idx = i - start;
            int gapAdjIdx = gapAdjustedIdx[idx];
            if (gapAdjIdx != (expectedGapAdjustedIdx)) {
                gaps.add(new Gap(i, expectedGapAdjustedIdx, gapAdjIdx));
            }
            expectedGapAdjustedIdx = gapAdjIdx + 1;
        }


        refSeq = new MASequence(bases.get(refSeqId));
        //assert ((refSeq != null) && (refSeq.length() > 0));
        //computeCoordinates();

        this.alignedSequences = new LinkedHashMap();
        for (Map.Entry<String, String> entry : bases.entrySet()) {
            String spId = entry.getKey();
            String b = entry.getValue();
            alignedSequences.put(spId, new MASequence(b));
        }
    }

    /**
     * @return the gaps
     */
    public List<Gap> getGaps() {
        return gaps;
    }

    /**
     * Class description
     *
     * @author Enter your name here...
     * @version Enter version here..., 09/01/22
     */
    public static class Gap {

        int position;
        int startIdx;
        int endIdx;

        /**
         * Constructs ...
         *
         * @param position
         * @param startIdx
         */
        public Gap(int position, int startIdx, int endIdx) {
            this.position = position;
            this.startIdx = startIdx;
            this.endIdx = endIdx;
        }

        public int getPosition() {
            return position;
        }

        /**
         * Method description
         *
         * @return
         */
        public int getCount() {
            return endIdx - startIdx;
        }

        /**
         * Method description
         *
         * @return
         */
        public String toString() {
            return position + "(" + (endIdx - startIdx + 1) + ")";
        }
    }

    public class MASequence {

        String bases;

        public MASequence(String bases) {
            this.bases = bases;
        }

        public char getGapAdjustedBase(int refCoord) {

            int relCoord = refCoord - start;
            int idx = gapAdjustedIdx[relCoord];
            return (idx > 0 && idx < bases.length()) ? bases.charAt(idx) : (char) 0;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }
    }
}
