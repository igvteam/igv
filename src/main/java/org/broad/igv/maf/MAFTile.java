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
