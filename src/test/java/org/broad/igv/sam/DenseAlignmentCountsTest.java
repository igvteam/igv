package org.broad.igv.sam;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.TestUtils;
import org.junit.Test;
import java.util.Iterator;
import static org.junit.Assert.*;

/**
 * @author jacob
 * @date 2013-Oct-31
 */
public class DenseAlignmentCountsTest extends AbstractHeadlessTest {

    /**
     * Test for IGV-2047, make sure we guard against the boundaries properly
     * @throws Exception
     */
    @Test
    public void testGetMaxCount_smallInterval() throws Exception {
        int fullIntervals = 0;
        int extraLength = DenseAlignmentCounts.MAX_COUNT_INTERVAL / 4 + 2;
        tstGetMaxCount(0, fullIntervals, extraLength);
    }

    /**
     * More normal case of largish interval
     * @throws Exception
     */
    @Test
    public void testGetMaxCount_LargeInterval() throws Exception {
        int fullIntervals = 10;
        int extraLength = DenseAlignmentCounts.MAX_COUNT_INTERVAL / 4 + 2;
        tstGetMaxCount(0, fullIntervals, extraLength);
    }

    @Test
    public void testEqBase() throws Exception {
        String chr = "chr1";
        int start = 59300;
        int end = 59400;
        String path = TestUtils.DATA_DIR + "bam/squeeze.sam";
        DenseAlignmentCounts daCounts = new DenseAlignmentCounts(start, end, null);

        AlignmentReader reader = AlignmentReaderFactory.getReader(path, false);
        Iterator<Alignment> iter = reader.iterator();
        while (iter.hasNext()) {
            Alignment al = iter.next();
            daCounts.incCounts(al);
        }
        int count = daCounts.getCount(59303, (byte) '=');
        assertEquals(3, count);

    }

    private void tstGetMaxCount(int start, int fullIntervals, int extraLength) {
        int mci = DenseAlignmentCounts.MAX_COUNT_INTERVAL;
        int end = start + fullIntervals * mci + extraLength;
        DenseAlignmentCounts daCounts = new DenseAlignmentCounts(start, end, null);

        int[] queryStarts = {start, start + 50, start + mci, start + mci + 50, start, start, start, start + 10};
        int[] queryLengths = {50, 50, 50, 50, (int) (mci * 1.5), 2 * mci, (int) (2.5 * mci), end - start + 10};
        for (int ii = 0; ii < queryStarts.length; ii++) {
            int queryStart = queryStarts[ii];
            int queryEnd = queryStart + queryLengths[ii];
            daCounts.getMaxCount(start, queryEnd);
        }

    }
}
