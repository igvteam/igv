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

package org.broad.igv.server;

import org.broad.igv.util.HttpUtils;
import org.junit.AfterClass;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.ProtocolException;
import java.net.URL;


public class ByteRangeTest {


    String urlString = "http://www.broadinstitute.org/igvdata/annotations/seq/hg19/chr1.txt";

    @Test
    public void testGetContentLength() throws IOException {
        // Open an input stream just to check permissions
        HttpURLConnection conn = null;
        try {
            conn = (HttpURLConnection) (new URL(urlString)).openConnection();
            String contentLength =  conn.getHeaderField("Content-length");
            assertEquals("249250621", contentLength);

        }
        finally {

            if (conn != null) {
                conn.disconnect();  // <- this really doesn't do anything (see Sun documentation)
            }
        }
    }

    @Test
    public void test2() throws Exception {
        assertTrue(HttpUtils.testByteRange("www.broadinstitute.org"));
    }

    @Test
    public void test1() throws IOException {

        int[] expectedBytes = {65, 67, 84, 65, 65, 71, 67, 65, 67, 65, 67};

        String byteRange = "bytes=" + 100000 + "-" + 100010;

        HttpURLConnection conn = (HttpURLConnection) (new URL(urlString)).openConnection();
        conn.setRequestMethod("GET");
        conn.setRequestProperty("Connection", "close");
        conn.setRequestProperty("Range", byteRange);

        InputStream is = null;
        try {
            is = conn.getInputStream();
            BufferedInputStream bis = new BufferedInputStream(is);

            for (int i = 0; i < expectedBytes.length; i++) {
                assertEquals("Testing index " + i, expectedBytes[i], bis.read());
            }

            assertEquals(-1, bis.read());
        }
        finally {
            if (is != null) {
                is.close();
            }
        }
    }
}
