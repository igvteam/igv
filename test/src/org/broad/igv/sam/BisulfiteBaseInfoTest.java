package org.broad.igv.sam;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Strand;
import org.junit.Test;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

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
        AlignmentBlock block = new AlignmentBlock(0, read, null, testAlignment);
        BisulfiteBaseInfo bbi = new BisulfiteBaseInfo(ref, block, context);

    }



    static class TestAlignment extends AbstractAlignment {

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
        public String getCigarString() {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
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

        public Strand getFragmentStrand(int read) {
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public boolean isProperPair() {
            return false;  //To change body of implemented methods use File | Settings | File Templates.
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
        private void setPairStrands() {

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
    }
}
