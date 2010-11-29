/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

/**
 * Timing test for S3 storage files
 */
public class TestS3Access {

    public static void main(String[] args) {
        try {
            testAsciiDownload();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    static void testAsciiDownload() throws IOException {

        String broadURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/refGene.txt";
        String s3URL = "http://s3.amazonaws.com/genomespace-igv/refGene.txt";

        try {


            for (int i = 0; i < 10; i++) {
                for (String url : Arrays.asList(broadURL, s3URL)) {

                    BufferedReader reader = openBufferedReader(url);
                    int nLines = 0;
                    long t0 = System.currentTimeMillis();
                    while (reader.readLine() != null) {
                        nLines++;
                    }
                    System.out.println(url + " read " + nLines + " in " + (System.currentTimeMillis() - t0) + " ms ");
                    reader.close();
                }
            }
        }
        finally {
        }
    }

    public static BufferedReader openBufferedReader(String pathOrUrl) throws IOException {

        BufferedReader reader = null;

        if (IGVHttpUtils.isURL(pathOrUrl)) {
            URL url = new URL(pathOrUrl);
            URLConnection connection = IGVHttpUtils.openConnection(url);
            reader = new BufferedReader(new InputStreamReader(connection.getInputStream()));
        } else {
            File file = new File(pathOrUrl);

            FileInputStream fileInput = new FileInputStream(file);
            if (file.getName().endsWith("gz")) {
                GZIPInputStream in = new GZIPInputStream(fileInput);
                reader = new BufferedReader(new InputStreamReader(in));
            } else {
                reader = new BufferedReader(new InputStreamReader(fileInput));
            }
        }

        return reader;
    }

}
