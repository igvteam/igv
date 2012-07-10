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
