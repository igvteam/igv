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
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


/**
 * Test of class HttpUtils (TODO -- need more tests!!)
 */
public class HttpUtilsTest extends AbstractHeadlessTest{


    static String broadURLString = "http://www.broadinstitute.org/igvdata/annotations/seq/hg19/chr1.txt";
    static String genericURLString = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz";

    static String noRangeHeaderSupportString = "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM714693&format=file&file=GSM714693%5Fhg19%5FwgEncodeGisDnaPetK562F1kAln%2Ebam";

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
     * Test of the Range test for a non-Broad URL.  Note if running this test behind proxies that strip
     * range headers the assertion should be false.
     *
     * @throws Exception
     */
    @Test
    public void testGenericURL() throws Exception {
        final URL url = new URL(genericURLString);
        String acceptsRangesValue = HttpUtils.getInstance().getHeaderField(url, "Accept-Ranges");
        boolean acceptsRanges = acceptsRangesValue != null && acceptsRangesValue.contains("bytes");
        assertEquals(acceptsRanges, HttpUtils.getInstance().useByteRange(url));
    }

    /**
     * Test that when we include the Range header, {@link HttpUtils#openConnection(java.net.URL, java.util.Map, String) we detect the absence of the range response properly
     *
     * @throws Exception
     */
    @Test
    public void testFailedRangeRequest() throws Exception {
        int start = 1;
        int end = 100;
        int numBytes = end - start + 1;
        String byteRange = "bytes=" + start + "-" + end;
        Map<String, String> params = new HashMap();
        params.put("Range", byteRange);

        HttpURLConnection conn = HttpUtils.getInstance().openConnection(new URL(noRangeHeaderSupportString), params);
        boolean rangeRequestedNotReceived = HttpUtils.getInstance().isExpectedRangeMissing(conn, params);

        assertTrue(rangeRequestedNotReceived);

        InputStream input = conn.getInputStream();
        byte b;
        int count=0, cc;
        while( (cc = input.read()) >= 0){
            count += cc;
            if(count > numBytes) break;
        }

        assertTrue("Read fewer than expected bytes even though range header ignored", count > numBytes);
        input.close();
    }

    @Ignore
    @Test
    public void testReadByBlob() throws Exception{
        String serverURL = "serverURLhere";

//        String bamURL =  serverURL + "?file=CDNLtx0vjHKW19PP9E7A1zpmzSg1mrg/rw7ZkE+b8FPuoNLNNLNKyL/94CsspYExn2SAJ5Y0wgL9L2mwf/HsP4IwqsQbtAqVBLxtd2CVclTOkute7VHhzqFPJ2KCzWxUjM+Ecb5XfdGTpKAzya1dq/fAtJuIw8+NHQCPmMVbJWreQx1+6z3VKDT17Fy3RbXxL6X/CfQ/HlTcWFpQVe1p+5LgkojOVagWCImYNk/ErPzi8J2oEYPSm6ilOwDM6rGwHcO47qW8DncaPdf8ohpm/XZwAd+aAJwsqkBR689R+X175QCzmpOI07dHxuvzJ4HPlMwu2h2290QxVAJ8Ix5fVA==.bam";
//
//        InputStream is = HttpUtils.getInstance().openConnectionStream(new URL(bamURL));

        String bedblob = StringUtils.decodeURL("c44z7H5E1gDMSm49T7NyGix051qDS7AbgCqicZ%2FpFkLobpmCim95byvYICc5VT%2Bv8Z%2FzE2gHWZkboBuME9eLxjEsfiO4bwnqZGP9fwoWXooK1LC8e3R5%2F6B9KyP9X3PR102PIApQASPfQGnYHqpBLifFPUbeRMqN%2BSxYi3h7udQJ8pli2QPEappIiOVWQ77cjJ6c0h12me6dT81fVrYT1E5CGHpnfUarIWCADRySVUfxqrwADKpnpaozMiWebh4OaWr5QfWHuG%2F%2F%2FhwVs7YxJaAR9S6pMqqfk213JofydHJjUimkv2V8tM3RJk3Q2y7CZ3Dk8X8wLiAJLfTIaXyreQ%3D%3D");

        String bedURL = serverURL + "?file=" + bedblob + ".bed";
        InputStream is = HttpUtils.getInstance().openConnectionStream(new URL(bedURL));

        BufferedReader br = new BufferedReader(new InputStreamReader(is));

        String line = null;
        int count = 0;
        while((line = br.readLine()) != null){
            System.out.println(line);
            count++;
        }
        System.out.println(count + " lines read");

    }

}
