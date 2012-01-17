package org.broad.igv.util;

import org.junit.Test;

import java.io.File;
import java.net.URL;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertTrue;

/**
 * @author Jim Robinson
 * @date 1/17/12
 */
public class HttpUtilsTest {
    @Test
    public void testGetContentsAsString() throws Exception {

    }

    @Test
    public void testGetContentLength() throws Exception {

    }

    @Test
    public void testCompareResources() throws Exception {

        File localFile = new File("test/data/bed/test.bed");
        assertTrue(localFile.exists());

        // To make this test robust we must compare the local file time to the known remote one
        long remoteModifedTime = 1326815079000l;
        boolean expectedValue = localFile.lastModified() >= remoteModifedTime;

        URL remoteURL = new URL("http://www.broadinstitute.org/igvdata/test/bed/test.bed");
        boolean comp = HttpUtils.getInstance().compareResources(localFile, remoteURL);

        assertEquals(expectedValue, comp);
    }


}
