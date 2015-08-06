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

import org.junit.Ignore;
import org.junit.Test;

import java.io.*;
import java.net.URL;
import java.util.Arrays;

import static org.junit.Assert.assertTrue;

/**
 * Timing test for S3 storage files
 */
@Ignore //Not sure about URLs
public class TestS3Access {

    public static void main(String[] args) {
        try {
            testAsciiDownloadStatic();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Test
    public void testAsciiDownload() throws Exception{
        testAsciiDownloadStatic();
    }
    public static void testAsciiDownloadStatic() throws IOException {

        String broadURL = "http://www.broadinstitute.org/igvdata/annotations/hg18/refGene.txt";
        String s3URL = "http://s3.amazonaws.com/genomespace-igv/refGene.txt";

        try {


            for (int i = 0; i < 10; i++) {
                for (String url : Arrays.asList(broadURL, s3URL)) {

                    InputStream inputStream = HttpUtils.getInstance().openConnectionStream(new URL(url));
                    BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
                    int nLines = 0;
                    long t0 = System.currentTimeMillis();
                    while (reader.readLine() != null) {
                        nLines++;
                    }
                    //System.out.println(url + " read " + nLines + " in " + (System.currentTimeMillis() - t0) + " ms ");
                    reader.close();
                    assertTrue(nLines > 0);
                }
            }
        } finally {
        }
    }
}
