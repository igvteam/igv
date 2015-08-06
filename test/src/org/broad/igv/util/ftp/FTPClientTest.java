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

import htsjdk.samtools.util.ftp.FTPClient;
import htsjdk.samtools.util.ftp.FTPReply;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 * @date Oct 30, 2010
 */
public class FTPClientTest {

    static String host = "ftp.broadinstitute.org";
    static String file = "/pub/igv/TEST/test.txt";
    static byte[] expectedBytes = "abcdefghijklmnopqrstuvwxyz\n".getBytes();
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
    public void testLogin() throws Exception {
        FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
        assertTrue(reply.isSuccess());

    }

    @Test
    public void testPasv() throws Exception {

        try {
            FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
            assertTrue(reply.isSuccess());

            reply = client.pasv();
            assertTrue(reply.isSuccess());
        }

        finally {
            client.closeDataStream();
        }
    }

    @Test
    public void testDownload() throws Exception {
        try {
            FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
            assertTrue(reply.isSuccess());

            reply = client.binary();
            assertTrue(reply.isSuccess());

            reply = client.pasv();
            assertTrue(reply.isSuccess());

            reply = client.retr(file);
            assertTrue(reply.getCode() == 150);

            InputStream is = client.getDataStream();
            int idx = 0;
            int b;
            while ((b = is.read()) >= 0) {
                assertEquals(expectedBytes[idx], (byte) b);
                idx++;
            }

        }

        finally {
            client.closeDataStream();
            FTPReply reply = client.retr(file);
            System.out.println(reply.getCode());
            assertTrue(reply.isSuccess());

        }
    }

    @Test
    public void testRest() throws Exception {
        try {
            FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
            assertTrue(reply.isSuccess());

            reply = client.binary();
            assertTrue(reply.isSuccess());

            reply = client.pasv();
            assertTrue(reply.isSuccess());

            final int restPosition = 5;
            client.setRestPosition(restPosition);

            reply = client.retr(file);
            assertTrue(reply.getCode() == 150);

            InputStream is = client.getDataStream();
            int idx = restPosition;
            int b;
            while ((b = is.read()) >= 0) {
                assertEquals(expectedBytes[idx], (byte) b);
                idx++;
            }

        }

        finally {
            client.closeDataStream();
            FTPReply reply = client.retr(file);
            System.out.println(reply.getCode());
            assertTrue(reply.isSuccess());

        }
    }

    @Test
    public void testMultiplePasv() throws Exception {

        try {
            FTPReply reply = client.login("anonymous", "igv@broadinstitute.org");
            assertTrue(reply.isSuccess());

            reply = client.pasv();
            assertTrue(reply.isSuccess());
            client.closeDataStream();

            reply = client.pasv();
            assertTrue(reply.isSuccess());
            client.closeDataStream();


        }

        finally {

        }
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

    private void restRetr(int restPosition, int length) throws IOException {

        try {

            if (client.getDataStream() == null) {
                FTPReply reply = client.pasv();
                assertTrue(reply.isSuccess());
            }

            client.setRestPosition(restPosition);

            FTPReply reply = client.retr(file);
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
            System.out.println(reply.getReplyString());


        }
    }
}
