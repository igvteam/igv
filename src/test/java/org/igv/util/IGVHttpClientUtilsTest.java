package org.igv.util;

import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date Jul 27, 2011
 */

public class IGVHttpClientUtilsTest {

    String url = "https://s3.amazonaws.com/igv.org.test/csi_test/chr10p.bam.csi";
    int byteCount = 105063;

    @Test
    public void testGetContentLength() throws IOException {

        assertEquals(byteCount, HttpUtils.getInstance().getContentLength(HttpUtils.createURL(url)));
    }

    @Test
    public void testExists() throws IOException {

        assertTrue("Resource unexpectedly does not exist", HttpUtils.getInstance().resourceAvailable(url));

        url = "http://nosuchserver/genomes/hg18.genome";
        assertFalse("Resource unexpectedly found", HttpUtils.getInstance().resourceAvailable(url));

        url = "http://igvdata.broadinstitute.org/nosuchfile.txt";
        assertFalse(HttpUtils.getInstance().resourceAvailable(url));
    }

}
