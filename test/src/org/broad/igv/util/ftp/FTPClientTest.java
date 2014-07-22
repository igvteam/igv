/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
