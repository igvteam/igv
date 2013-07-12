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
