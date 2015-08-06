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
