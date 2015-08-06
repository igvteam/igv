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
