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

package org.broad.igv.util;

import org.broad.igv.AbstractHeadlessTest;
import org.junit.Test;

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

import static org.junit.Assert.*;


/**
 * Test of class HttpUtils (TODO -- need more tests!!)
 */
public class HttpUtilsTest extends AbstractHeadlessTest {


    static String broadURLString = "https://igv-genepattern-org.s3.amazonaws.com/test/fasta/chr22.fa";

    @Test
    public void testSignedURLMatch() throws Exception {
        String aws = "https://amazonaws.com?X-Amz-Signature=foo";  //X-Amz-Signature"
        assertTrue(HttpUtils.isSignedURL(aws));
        String google = "https://google.com?X-Goog-Signature=bar";
        assertTrue(HttpUtils.isSignedURL(google));
        String noMatch = "https://www.google.com";
        assertFalse(HttpUtils.isSignedURL(noMatch));
    }


    @Test
    public void testGetContentLength() throws IOException {
        // Open an input stream just to check permissions
        HttpURLConnection conn = null;
        try {
            conn = (HttpURLConnection) (HttpUtils.createURL(broadURLString)).openConnection();
            String contentLength = conn.getHeaderField("Content-length");
            assertEquals("52330665", contentLength);

        } finally {

            if (conn != null) {
                conn.disconnect();  // <- this really doesn't do anything (see Sun documentation)
            }
        }
    }


    @Test
    public void testCacheControl() throws Exception {
        String headerValue = "no-cache, max-age=100";
        HttpUtils.CacheControl cc = HttpUtils.CacheControl.valueOf(headerValue);
        assertTrue(cc.isNoCache());
        assertEquals(100, cc.getMaxAge());
    }

    @Test
    public void testAccessTokenCache() throws MalformedURLException {

        try {
            // Exact match
            HttpUtils.getInstance().setAccessToken("foo", "bar.foo.com");
            String token = HttpUtils.getInstance().getCachedTokenFor(new URL("https://bar.foo.com/path"));
            assertEquals("foo", token);
            HttpUtils.getInstance().clearAccessTokens();

            // Wildcard match
            HttpUtils.getInstance().setAccessToken("foo", "*.foo.com");
            token = HttpUtils.getInstance().getCachedTokenFor(new URL("https://bar.foo.com/path"));
            assertEquals("foo", token);

            // Superceding match
            HttpUtils.getInstance().setAccessToken("foo2", "*.foo.com");
            token = HttpUtils.getInstance().getCachedTokenFor(new URL("https://bar.foo.com/path"));
            assertEquals("foo2", token);


            // Clear token
            HttpUtils.getInstance().clearAccessTokens();
            token = HttpUtils.getInstance().getCachedTokenFor(new URL("https://bar.foo.com/path"));
            assertNull(token);
            HttpUtils.getInstance().clearAccessTokens();

            // Match all hosts
            HttpUtils.getInstance().setAccessToken("foo", "");
            token = HttpUtils.getInstance().getCachedTokenFor(new URL("https://igv.org/path"));
            assertEquals("foo", token);
        } finally {
            HttpUtils.getInstance().clearAccessTokens();
        }
    }
}
