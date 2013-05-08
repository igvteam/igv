package org.broad.igv.track;

import org.broad.igv.Globals;
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
        TribbleFeatureSource featuresSource = new TribbleFeatureSource(file, null, false);
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
        TribbleFeatureSource featuresSource = new TribbleFeatureSource(file, null, false);

        for (int tri = 0; tri < trials; tri++) {
            featuresSource.getFeatures("chr2", 178709699, 178711955);
            int curFileHandles = TestUtils.getNumberOpenFileHandles();
            String msg = "Number of open file handles deviates too much from base:\n";
            msg += "base: " + baseFileHandles + " current: " + curFileHandles + " trial number: " + tri;
            Assert.assertTrue(msg, curFileHandles - baseFileHandles <= maxDiff);
        }
    }
}
