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

package org.broad.igv.track;

import org.broad.igv.Globals;
import org.broad.igv.util.ResourceLocator;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;
import org.broad.igv.util.TestUtils;


/**
 * @author jrobinso
 *         Date: 5/8/13
 *         Time: 8:32 AM
 */
public class TribbleFeatureSourceTest {


    /**
     * Test release of file handles.  In previous versions of IGV file handles were not closed after a
     * tribble indexed query, leading to a "too many file handles" error.
     *
     * @throws Exception
     */
    @Test
    public void testCloseFileHandles() throws Exception {
        String file = org.broad.igv.util.TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";
        TribbleFeatureSource featuresSource = TribbleFeatureSource.getFeatureSource(new ResourceLocator(file), null, false);
        for (int i = 0; i < 10000; i++) {
            featuresSource.getFeatures("chr2", 178709699, 178711955);
        }
    }



    @Ignore   // Ignored, getting Cannot run program "lsof" error on test server. The test is redundant to testCloseFileHandles()
    @Test
    public void testFileHandleNumberNonincreasing() throws Exception {
        int trials = 100;
        int maxDiff = 2;
        Assume.assumeTrue(Globals.IS_LINUX || Globals.IS_MAC);
        int baseFileHandles = TestUtils.getNumberOpenFileHandles();

        String file = org.broad.igv.util.TestUtils.DATA_DIR + "bed/Unigene.sample.sorted.bed";
        TribbleFeatureSource featuresSource = TribbleFeatureSource.getFeatureSource(new ResourceLocator(file), null, false);

        for (int tri = 0; tri < trials; tri++) {
            featuresSource.getFeatures("chr2", 178709699, 178711955);
            int curFileHandles = TestUtils.getNumberOpenFileHandles();
            String msg = "Number of open file handles deviates too much from base:\n";
            msg += "base: " + baseFileHandles + " current: " + curFileHandles + " trial number: " + tri;
            Assert.assertTrue(msg, curFileHandles - baseFileHandles <= maxDiff);
        }
    }
}
