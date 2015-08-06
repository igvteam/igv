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

package org.broad.igv.util.ftp;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.stream.IGVSeekableFTPStream;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;

/**
 * User: jrobinso
 * Date: Apr 13, 2010
 */
public class SeekableFTPStreamTest extends AbstractHeadlessTest {

    static String url = "ftp://ftp.broadinstitute.org/distribution/igv/TEST/cpgIslands.hg18.bed";

    @Test
    public void testRead() throws IOException {

        byte[] buffer = new byte[100];

        IGVSeekableFTPStream stream = new IGVSeekableFTPStream(new URL(url));

        for (int j = 0; j < 10; j++) {
            stream.seek(20 + j * 100);

            stream.read(buffer, 0, buffer.length);

            //for (int i = 0; i < 10; i++) {
            //System.out.println(buffer[i]);
            //}
        }
        stream.close();
    }

    @Test
    public void testRead2() throws IOException {

        byte[] buffer = new byte[100];

        URL u = new URL(url);
        URLConnection conn = u.openConnection();

        conn.connect();

        //PrintWriter pw = new PrintWriter(new OutputStreamWriter(conn.getOutputStream()));
        //pw.close();

        InputStream stream = conn.getInputStream();
        for (int j = 0; j < 1; j++) {

            stream.read(buffer, 0, buffer.length);

            //for (int i = 0; i < 10; i++) {
            //System.out.println(buffer[i]);
            //}
        }

        stream.close();


    }
}
