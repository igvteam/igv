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


}
