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

package org.broad.igv.util.ftp;

import junit.framework.TestCase;
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
public class SeekableFTPStreamTest extends TestCase {

    static String url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase1/data/HG00099/alignment/HG00099.mapped.SOLID.bfast.GBR.low_coverage.20101123.bam.bai";

    @Test
    public void testRead() throws IOException {

        byte[] buffer = new byte[100];

        System.out.println((new URL(url)).getHost());
        IGVSeekableFTPStream stream = new IGVSeekableFTPStream(new URL(url));

        for (int j = 0; j < 10; j++) {
            stream.seek(20 + j * 100);

            stream.read(buffer, 0, buffer.length);

            for (int i = 0; i < 10; i++) {
                //System.out.println(buffer[i]);
            }
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

            for (int i = 0; i < 10; i++) {
                //System.out.println(buffer[i]);
            }
        }

        stream.close();


    }
}
