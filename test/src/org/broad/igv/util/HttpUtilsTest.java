/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.util;

import org.broad.igv.util.HttpUtils;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;


/**
 * Test of class HttpUtils (TODO -- need more tests!!)
 */
public class HttpUtilsTest {


    String broadURLString = "http://www.broadinstitute.org/igvdata/annotations/seq/hg19/chr1.txt";
    String genericURLString = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz";

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
