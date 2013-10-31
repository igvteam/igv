/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.sam;

import org.broad.igv.AbstractHeadlessTest;
import org.junit.Test;

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

    private void tstGetMaxCount(int start, int fullIntervals, int extraLength){
        int mci = DenseAlignmentCounts.MAX_COUNT_INTERVAL;
        int end = start + fullIntervals*mci + extraLength;
        DenseAlignmentCounts daCounts = new DenseAlignmentCounts(start, end, null);

        int[] queryStarts = {start, start + 50, start + mci, start + mci + 50, start,           start, start, start + 10};
        int[] queryLengths ={50,    50,         50,          50,               (int) (mci*1.5), 2*mci, (int)(2.5*mci), end-start + 10};
        for(int ii=0; ii < queryStarts.length; ii++){
            int queryStart = queryStarts[ii];
            int queryEnd = queryStart + queryLengths[ii];
            daCounts.getMaxCount(start, queryEnd);
        }

    }
}
