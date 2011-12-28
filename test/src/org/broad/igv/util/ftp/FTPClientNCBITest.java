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


import org.broad.tribble.util.ftp.FTPClient;
import org.broad.tribble.util.ftp.FTPReply;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.SocketException;
import java.net.SocketTimeoutException;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 * @date Oct 30, 2010
 */
public class FTPClientNCBITest {

    static String host = "ftp-trace.ncbi.nih.gov";

    private static String dirs = "1000genomes/ftp/data/HG00099/alignment/";
    private static String filename = "HG00099.chrom20.SOLID.bfast.GBR.low_coverage.20101123.bam.bai";
    private static String file = dirs + "/" + filename;
    private static byte[] expectedBytes;
    FTPClient client;

    private void loadLocalFile() throws IOException {
        String path = "test/data/bam/" + filename;
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(path));
        expectedBytes = new byte[bis.available()];
        bis.read(expectedBytes);
    }

    @Before
    public void setUp() throws IOException {
        client = new FTPClient();
        FTPReply reply = client.connect(host);
        assertTrue(reply.isSuccess());

        loadLocalFile();
    }

    @After
    public void tearDown() {
        //System.out.println("Disconnecting");
        client.disconnect();
        expectedBytes = null;
    }


    @Test
    public void testRest() throws Exception {
        FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
        assertTrue(reply.isSuccess());

        reply = client.binary();
        assertTrue(reply.isSuccess());

        restRetr(10, 10);
    }


    @Test
    public void testMultipleRest() throws Exception {
        FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
        assertTrue(reply.isSuccess());

        reply = client.binary();
        assertTrue(reply.isSuccess());

        restRetr(5, 10);
        restRetr(2, 10);
        restRetr(15, 10);

    }

    private void restRetr(int restPosition, int length) throws IOException, InterruptedException {
        boolean getReply = true;
        try {

            FTPReply reply = client.pasv();
            assertTrue(reply.isSuccess());

            client.setRestPosition(restPosition);

            client.retr(file);
            //assertTrue(reply.getCode() == 150);

            InputStream is = client.getDataStream();

            byte[] buffer = new byte[length];
            is.read(buffer);


            for (int i = 0; i < length; i++) {
                System.out.print((char) buffer[i]);
                assertEquals(expectedBytes[i + restPosition], buffer[i]);
            }
        } catch (SocketException e) {
            //Error contacting server
            //This isn't strictly necessary, but avoids trying
            //to contact server in finally clause
            getReply = false;
        } finally {
            client.closeDataStream();
            if (getReply) {
                //If we were able
                client.getReply();
            }
        }

        if (!getReply) {
            throw new SocketTimeoutException("Could not contact host");
        }
    }
}

