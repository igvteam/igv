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

import org.broad.igv.PreferenceManager;
import org.broad.igv.exceptions.HttpResponseException;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.net.URL;

import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date Jul 27, 2011
 */

public class IGVHttpClientUtilsTest {

    private static String hg18URL = "http://www.broadinstitute.org/igvdata/test/hg18.unittest.genome";
    private static int hg18bytes = 3617644;

    @Test
    public void testGetContentLength() throws IOException {

        String url = hg18URL;
        assertEquals(hg18bytes, HttpUtils.getInstance().getContentLength(new URL(url)));
    }

    @Test
    public void testExists() throws IOException {
        URL url = new URL(hg18URL);
        assertTrue("Resource unexpectedly does not exist", HttpUtils.getInstance().resourceAvailable(url));

        url = new URL("http://nosuchserver/genomes/hg18.genome");
        assertFalse("Resource unexpectedly found", HttpUtils.getInstance().resourceAvailable(url));

        url = new URL("http://igvdata.broadinstitute.org/nosuchfile.txt");
        assertFalse(HttpUtils.getInstance().resourceAvailable(url));
    }

    /**
     * This test will only work if on the Broad intranet, as that's where the test proxy server is.
     *
     * @throws IOException
     */
    @Ignore
    @Test
    public void testProxy() throws IOException {

        final URL testURL = new URL(hg18URL);

        PreferenceManager mgr = PreferenceManager.getInstance();
        mgr.override(PreferenceManager.PROXY_HOST, "igvdev01.broadinstitute.org");
        mgr.override(PreferenceManager.PROXY_PORT, "3128");
        mgr.override(PreferenceManager.PROXY_USER, "proxytest");
        String enc_pword = Utilities.base64Encode("test@123");
        mgr.override(PreferenceManager.PROXY_PW, enc_pword);
        mgr.override(PreferenceManager.USE_PROXY, "true");
        mgr.override(PreferenceManager.PROXY_AUTHENTICATE, "true");
        HttpUtils.getInstance().updateProxySettings();

        long contentLength = 0;
        try {
            contentLength = HttpUtils.getInstance().getContentLength(testURL);
            assertEquals(hg18bytes, contentLength);
        } catch (IOException e) {
            System.out.println("Proxy unreachable.  Skipping proxy test");
            return;
        }

        // Now try to get a file not on the squid "allowed" domains to verify requests are going through the proxy
        // This should fail and return -1 for the content length
        try {
            contentLength = HttpUtils.getInstance().getContentLength(new URL("http://www.boston.com"));
            junit.framework.Assert.fail("Proxy test is apparently bypassing proxy");
        } catch (HttpResponseException e) {
            // This is expected
            assertEquals(403, e.getStatusCode());
        }


    }
}
