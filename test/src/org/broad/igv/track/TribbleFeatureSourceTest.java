package org.broad.igv.track;

import org.junit.Test;

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
        String file =  org.broad.igv.util.TestUtils.DATA_DIR +  "bed/Unigene.noheader.sorted.bed";
        TribbleFeatureSource featuresSource = new TribbleFeatureSource(file, null, false);
        for(int i=0; i<100000; i++) {
            featuresSource.getFeatures("chr2", 178709699, 178711955);
        }
    }
}
