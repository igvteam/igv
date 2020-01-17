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
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static spark.Spark.get;


/**
 * Test of class HttpUtils (TODO -- need more tests!!)
 */
public class HttpUtilsTest extends AbstractHeadlessTest {


    static String broadURLString = "http://data.broadinstitute.org/igvdata/annotations/seq/hg19/chr1.txt";
    static String genericURLString = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz";

    static String noRangeHeaderSupportString = "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM714693&format=file&file=GSM714693%5Fhg19%5FwgEncodeGisDnaPetK562F1kAln%2Ebam";

    @Test
    public void testGetContentLength() throws IOException {
        // Open an input stream just to check permissions
        HttpURLConnection conn = null;
        try {
            conn = (HttpURLConnection) (HttpUtils.createURL(broadURLString)).openConnection();
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
        assertTrue(HttpUtils.getInstance().useByteRange(HttpUtils.createURL(broadURLString)));
    }


    /**
     * Test of the Range test for a non-Broad URL.  Note if running this test behind proxies that strip
     * range headers the assertion should be false.
     *
     * @throws Exception
     */
    @Test
    public void testGenericURL() throws Exception {
        final URL url = HttpUtils.createURL(genericURLString);
        String acceptsRangesValue = HttpUtils.getInstance().getHeaderField(url, "Accept-Ranges");
        boolean acceptsRanges = acceptsRangesValue != null && acceptsRangesValue.contains("bytes");
        assertEquals(acceptsRanges, HttpUtils.getInstance().useByteRange(url));
    }


    @Ignore
    @Test
    public void testReadByBlob() throws Exception {
        String serverURL = "serverURLhere";

//        String bamURL =  serverURL + "?file=CDNLtx0vjHKW19PP9E7A1zpmzSg1mrg/rw7ZkE+b8FPuoNLNNLNKyL/94CsspYExn2SAJ5Y0wgL9L2mwf/HsP4IwqsQbtAqVBLxtd2CVclTOkute7VHhzqFPJ2KCzWxUjM+Ecb5XfdGTpKAzya1dq/fAtJuIw8+NHQCPmMVbJWreQx1+6z3VKDT17Fy3RbXxL6X/CfQ/HlTcWFpQVe1p+5LgkojOVagWCImYNk/ErPzi8J2oEYPSm6ilOwDM6rGwHcO47qW8DncaPdf8ohpm/XZwAd+aAJwsqkBR689R+X175QCzmpOI07dHxuvzJ4HPlMwu2h2290QxVAJ8Ix5fVA==.bam";
//
//        InputStream is = HttpUtils.getInstance().openConnectionStream(HttpUtils.createURL(bamURL));

        String bedblob = StringUtils.decodeURL("c44z7H5E1gDMSm49T7NyGix051qDS7AbgCqicZ%2FpFkLobpmCim95byvYICc5VT%2Bv8Z%2FzE2gHWZkboBuME9eLxjEsfiO4bwnqZGP9fwoWXooK1LC8e3R5%2F6B9KyP9X3PR102PIApQASPfQGnYHqpBLifFPUbeRMqN%2BSxYi3h7udQJ8pli2QPEappIiOVWQ77cjJ6c0h12me6dT81fVrYT1E5CGHpnfUarIWCADRySVUfxqrwADKpnpaozMiWebh4OaWr5QfWHuG%2F%2F%2FhwVs7YxJaAR9S6pMqqfk213JofydHJjUimkv2V8tM3RJk3Q2y7CZ3Dk8X8wLiAJLfTIaXyreQ%3D%3D");

        String bedURL = serverURL + "?file=" + bedblob + ".bed";
        InputStream is = HttpUtils.getInstance().openConnectionStream(HttpUtils.createURL(bedURL));

        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        String line = null;
        int count = 0;
        while ((line = br.readLine()) != null) {
            System.out.println(line);
            count++;
        }
        System.out.println(count + " lines read");

    }

    @Test
    public void testCacheControl() throws Exception {
        String headerValue = "no-cache, max-age=100";
        HttpUtils.CacheControl cc = HttpUtils.CacheControl.valueOf(headerValue);
        assertTrue(cc.isNoCache());
        assertEquals(100, cc.getMaxAge());
    }

    public class RunnableSparkHttp implements Runnable {

        /// counters incremented on each request:

        // permanent redirect src
        public int permSrcCt = 0;
        // permanent redirect dest
        public int permDestCt = 0;
        // temporary redirect src
        public int tempSrcCt = 0;
        // temporary redirect dest
        public int tempDestCt = 0;

        public void run() {
            System.out.println("run thing");

            get("/perm_redir_src", (req, res) ->
            {
                // increment a counter for tests to inspect
                this.permSrcCt += 1;
                res.status(301); // permanent redirect
                res.header("Location", "http://localhost:4567/perm_redir_dest");
                return "redirecting";
            });

            get("/perm_redir_dest", (req, res) -> {
                this.permDestCt += 1;
                return "done";
            });

            get("/temp_redir_src", (req, res) ->
            {
                this.tempSrcCt += 1;
                res.status(302); // temporary redirect
                res.header("Location", "http://localhost:4567/temp_redir_dest");
                res.header("Cache-Control", "max-age=1"); // expire in 1 second
                return "redirecting";
            });

            get("/temp_redir_dest", (req, res) -> {
                this.tempDestCt += 1;
                return "done";
            });

        }

    }

    @Test
    public void testRedirectCache() throws Exception {

        RunnableSparkHttp server = new RunnableSparkHttp();
        Thread thread = new Thread(server);
        thread.run();

        // give the test server a moment to start up
        Thread.sleep(500);

        HttpURLConnection conn = null;

        // test a URL that redirects permanently
        assertEquals(server.permSrcCt, 0);
        assertEquals(server.permDestCt, 0);
        conn = HttpUtils.getInstance().openConnection(new URL("http://localhost:4567/perm_redir_src"), null);
        assertEquals(conn.getResponseCode(), 200);
        assertEquals(server.permSrcCt, 1);
        assertEquals(server.permDestCt, 1);
        assertEquals(HttpUtils.getInstance().openConnection(new URL("http://localhost:4567/perm_redir_src"), null)
                .getResponseCode(), 200);
        // because redirect was cached, the source wasn't requested again, but the destination was
        assertEquals(server.permSrcCt, 1);
        assertEquals(server.permDestCt, 2);


        // now test a URL that redirects but quickly expires
        assertEquals(server.tempSrcCt, 0);
        assertEquals(server.tempDestCt, 0);
        conn = HttpUtils.getInstance().openConnection(new URL("http://localhost:4567/temp_redir_src"), null);
        assertEquals(conn.getResponseCode(), 200);
        assertEquals(server.tempSrcCt, 1);
        assertEquals(server.tempDestCt, 1);
        // sleep for 2s to ensure that cache has expired
        Thread.sleep(2_000);
        assertEquals(HttpUtils.getInstance().openConnection(new URL("http://localhost:4567/temp_redir_src"), null)
                .getResponseCode(), 200);
        // because redirect was cached, the source wasn't requested again, but the destination was
        assertEquals(server.tempSrcCt, 2);
        assertEquals(server.tempDestCt, 2);

    }
}
