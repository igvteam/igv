/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.util;

import org.junit.Test;

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


/**
 * Test of class HttpUtils (TODO -- need more tests!!)
 */
public class HttpUtilsTest {


    static String broadURLString = "http://www.broadinstitute.org/igvdata/annotations/seq/hg19/chr1.txt";
    static String genericURLString = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz";

    @Test
    public void testGetContentLength() throws IOException {
        // Open an input stream just to check permissions
        HttpURLConnection conn = null;
        try {
            conn = (HttpURLConnection) (new URL(broadURLString)).openConnection();
            String contentLength = conn.getHeaderField("Content-length");
            assertEquals("249250621", contentLength);

        } finally {

            if (conn != null) {
                conn.disconnect();  // <- this really doesn't do anything (see Sun documentation)
            }
        }
    }

    /**
     * Test of the "byte range" test for a Broad URL.  Note if running this test behind proxies that strip
     * range headers the assertion should be "false".
     *
     * @throws Exception
     */
    @Test
    public void testBroadURL() throws Exception {
        assertTrue(HttpUtils.getInstance().useByteRange(new URL(broadURLString)));
    }


    /**
     * Test of the "byte range" test for a non-Broad URL.  Note if running this test behind proxies that strip
     * range headers the assertion should be "false".
     *
     * @throws Exception
     */    @Test
    public void testGenericURL() throws Exception {
        final URL url = new URL(genericURLString);
        String acceptsRangesValue = HttpUtils.getInstance().getHeaderField(url, "Accept-Ranges");
        boolean acceptsRanges = acceptsRangesValue != null && acceptsRangesValue.contains("bytes");
        assertEquals(acceptsRanges, HttpUtils.getInstance().useByteRange(url));
    }

}
