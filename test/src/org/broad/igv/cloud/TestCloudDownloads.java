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

package org.broad.igv.cloud;


import org.broad.igv.util.HttpUtils;
import org.junit.Test;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * User: jrobinso
 * Date: Feb 4, 2010
 */
public class TestCloudDownloads {

    public TestCloudDownloads() {
    }

    @Test
    public void testSequenceByteRange() throws IOException {

        String broadURL = "http://www.broadinstitute.org/igvdata/annotations/seq/hg18/chr1.txt";
        String s3URL = "http://igv.broadinstitute.org/genomes/seq/hg18/chr1.txt";
        String cfURL = "http://igvdata.broadinstitute.org/genomes/seq/hg18/chr1.txt";

        int size = 100000;
        int start = Math.max(1, (int) (Math.random() * 2000000));
        for (String url : Arrays.asList(broadURL, s3URL, cfURL)) {
            final int nBytes = readByteRange(url, start, start + (size - 1));
            assertEquals(size, nBytes);
        }



    }


    public int readByteRange(String urlString, int start, int end) throws IOException {

        URL url = new URL(urlString);
        Map<String, String> reqProps = new HashMap();
        String byteRange = "bytes=" + start + "-" + end;
        reqProps.put("Range", byteRange);
        InputStream is = HttpUtils.getInstance().openConnectionStream(url, reqProps);
        BufferedInputStream bis = new BufferedInputStream(is);
        int b;
        int n = 0;

        while ((b = bis.read()) >= 0) {
            n++;
        }
        //System.out.println("Read " + n + " bytes");

        is.close();

        return n;
    }

}
