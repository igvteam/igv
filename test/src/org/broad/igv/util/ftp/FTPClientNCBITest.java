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

package org.broad.igv.util.ftp;

import org.broad.igv.util.ftp.FTPClient;
import org.broad.igv.util.ftp.FTPReply;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00099/alignment/HG00099.chrom1.SOLID.bfast.GBR.low_coverage.20100817.bam;
 * ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00099/alignment/HG00099.chrom1.SOLID.bfast.GBR.low_coverage.20100817.bam
 *
 * @author jrobinso
 * @date Oct 30, 2010
 */
public class FTPClientNCBITest {

    //static String host = "ftp.1000genomes.ebi.ac.uk";
    //static String file = "/vol1/ftp/data/HG00099/alignment/HG00099.chrom1.SOLID.bfast.GBR.low_coverage.20100817.bam";
    //ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00099/alignment/HG00099.chrom20.SOLID.bfast.GBR.low_coverage.20101123.bam
    static String host = "ftp-trace.ncbi.nih.gov";
    static String file = "/1000genomes/ftp/data/HG00099/alignment/HG00099.chrom1.SOLID.bfast.GBR.low_coverage.20100817.bam";
    static byte[] expectedBytes = {(byte) 31, (byte) 139, (byte) 8, (byte) 4, 0, 0, 0, 0, 0, (byte) 255, (byte) 6, 0,
            (byte) 66, (byte) 67, (byte) 2, 0, (byte) 240, (byte) 100, (byte) 204, (byte) 189, (byte) 123, (byte) 112, (byte) 99,
            (byte) 105, (byte) 118, (byte) 31, (byte) 134};
    FTPClient client;

    @Before
    public void setUp() throws IOException {
        client = new FTPClient();
        FTPReply reply = client.connect(host);
        assertTrue(reply.isSuccess());
    }

    @After
    public void tearDown() {
        System.out.println("Disconnecting");
        client.disconnect();
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

        try {

            FTPReply reply = client.pasv();
            assertTrue(reply.isSuccess());

            client.setRestPosition(restPosition);

            reply = client.retr(file);
            //assertTrue(reply.getCode() == 150);

            InputStream is = client.getDataStream();

            byte[] buffer = new byte[length];
            is.read(buffer);


            for (int i = 0; i < length; i++) {
                System.out.print((char) buffer[i]);
                assertEquals(expectedBytes[i + restPosition], buffer[i]);
            }
            System.out.println();
        }

        finally {

            client.closeDataStream();
            FTPReply reply = client.getReply();  // <== MUST READ THE REPLY

        }
    }
}

