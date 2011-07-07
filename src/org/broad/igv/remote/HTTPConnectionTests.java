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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.remote;

import org.broad.igv.util.IGVHttpClientUtils;
import org.broad.igv.util.IGVHttpUtils;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class HTTPConnectionTests {


    public static long getContentLength(URL url) throws IOException {

        return IGVHttpClientUtils.getContentLength(url);
    }

    public static void dumpHeaderFields() throws IOException {


   /*     URL url = new URL("http://www.broadinstitute.org/igv/resources/dataServerRegistry.txt");
        HttpURLConnection connection = (HttpURLConnection) IGVHttpUtils.openConnection(url);

        connection.setRequestMethod("HEAD");

        connection.connect();

        Map<String, List<String>> map = connection.getHeaderFields();

        for (Map.Entry<String, List<String>> entry : map.entrySet()) {
            System.out.print(entry.getKey() + ":\t");
            for (String v : entry.getValue()) {
                System.out.print(v + " ");
            }
            System.out.println();
        }


        connection.disconnect();     */
    }

    public static void getByteRange() throws IOException {


    /*    URL url = new URL("http://www.broadinstitute.org/igv/resources/dataServerRegistry.txt");
        //URL url = new URL("http://www.broadinstitute.org/~jrobinso/dataServerRegistry.txt");

        long len = getContentLength(url);

        HttpURLConnection connection = (HttpURLConnection) IGVHttpUtils.openConnection(url);
        connection.setRequestMethod("POST");

        String byteRange = "bytes=" + (len - 10) + "-" + (len - 1);
        System.out.println(byteRange);
        connection.setRequestProperty("Range", byteRange);

        Map<String, List<String>> map = connection.getHeaderFields();

        for (Map.Entry<String, List<String>> entry : map.entrySet()) {
            System.out.print(entry.getKey() + ":\t");
            for (String v : entry.getValue()) {
                System.out.print(v + " ");
            }
            System.out.println();
        }


        InputStream is = connection.getInputStream();
        char b;
        int n = 0;
        while ((b = (char) is.read()) >= 0 && n < 9) {
            System.out.print(b);
            n++;
        }

        connection.disconnect();   */
    }
}
