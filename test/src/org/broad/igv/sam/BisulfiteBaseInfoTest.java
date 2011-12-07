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
    }
}
