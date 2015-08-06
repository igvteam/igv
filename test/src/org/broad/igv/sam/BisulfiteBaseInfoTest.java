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

package org.broad.igv.sam;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.junit.Test;

/**
 * @author Jim Robinson
 * @date 12/6/11
 */
public class BisulfiteBaseInfoTest {

    @Test
    public void testContextIsMatching() throws Exception {

        AlignmentTrack.BisulfiteContext context =  AlignmentTrack.BisulfiteContext.CG;
        byte [] ref = {(byte) 'C', (byte) 'G'};
        byte [] read = {(byte) 'C', (byte) 'G'};
        TestAlignment testAlignment = new TestAlignment();
        AlignmentBlock block = new AlignmentBlock("", 0, read, null);
        BisulfiteBaseInfo bbi = new BisulfiteBaseInfo(ref, testAlignment, block, context);

    }



    static class TestAlignment extends SAMAlignment {

        boolean isNegativeStrand = false;
        private Strand firstOfPairStrand;
        private Strand secondOfPairStrand;

        public String getReadSequence() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getAlignmentStart() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public String getReadName() {
            return null;
        }

        @Override
        public int getMappingQuality() {
            return 0;
        }

        @Override
        public int getInferredInsertSize() {
            return 0;
        }

        @Override
        public String getCigarString() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public int getReadLength() {
            return 0;
        }

        @Override
        public boolean isMapped() {
            return false;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public boolean isPaired() {
            return false;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public boolean isFirstOfPair() {
            return false;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public boolean isSecondOfPair() {
            return false;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public boolean isDuplicate() {
            return false;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public int getAlignmentEnd() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public String getSample() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public Object getAttribute(String key) {
            return null;
        }

        public Strand getFragmentStrand(int read) {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public boolean isProperPair() {
            return false;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public boolean isSupplementary() {
            return false;
        }

        @Override
        public boolean isVendorFailedRead() {
            return false;
        }

        @Override
        public boolean isPrimary() {
            return false;
        }

        public void setStart(int start) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public void setEnd(int end) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        public LocusScore copy() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        /**
         * Return the start position
         */
        public int getStart() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        /**
         * Return the end position
         */
        public int getEnd() {
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        /**
         * Set the strands of first and second in pair.  Used for strand specific libraries to recover strand of
         * originating fragment.
         */
        protected void setPairStrands() {

            if (isPaired()) {
                if (isFirstOfPair()) {
                    firstOfPairStrand = getReadStrand();
                } else {
                    // If we have a mate, the mate must be the firstOfPair
                    ReadMate mate = getMate();
                    if (mate != null && mate.isMapped()) {
                        firstOfPairStrand = mate.getStrand();
                    } else {
                        // No Mate, or mate is not mapped, FOP strand is not defined
                        firstOfPairStrand = Strand.NONE;
                    }
                }

                if (isSecondOfPair()) {
                    secondOfPairStrand = isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
                } else {
                    ReadMate mate = getMate();
                    if (mate.isMapped() && isProperPair()) {
                        secondOfPairStrand = mate.isNegativeStrand() ? Strand.NEGATIVE : Strand.POSITIVE;
                    } else {
                        // No Mate, or mate is not mapped, FOP strand is not defined
                        secondOfPairStrand = Strand.NONE;
                    }
                }

            } else {
                // This alignment is not paired -- by definition "firstOfPair" is this alignment
                firstOfPairStrand = getReadStrand();
                secondOfPairStrand = Strand.NONE;
            }
        }

        @Override
        protected String getAttributeString(boolean truncate) {
            return "";
        }

        /**
         * Return the strand of the read marked "first-in-pair" for a paired alignment. This method can return
         * Strand.NONE if the end marked first is unmapped.
         *
         * @return strand of first-of-pair
         */
        public Strand getFirstOfPairStrand() {
           return firstOfPairStrand;
        }

        /**
         * Return the strand of the read marked "second-in-pair" for a paired alignment.  The strand is
         * undefined (Strand.NONE) for non-paired alignments
         *
         * @return strand of second-of-pair
         */
        public Strand getSecondOfPairStrand() {
            return secondOfPairStrand;
        }

        public boolean isNegativeStrand() {
            return false;
        }
    }
}
