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



import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;

import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 * @date Aug 31, 2010
 */
public class FTPUtils {

    public static boolean resourceAvailable(URL url) {
        InputStream is = null;
        try {
            URLConnection conn = url.openConnection();
            conn.setConnectTimeout(10000);
            conn.setReadTimeout(10000);
            is = conn.getInputStream();
            return (is.read() >= 0);

        } catch (IOException e) {
            return false;
        }
        finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        }
    }


    /**
     * Connect to an FTP server anonymously
     * 
     * @param host
     * @return
     * @throws IOException
     */
    public static FTPClient connect(String host)  throws IOException {

        FTPClient ftp = new FTPClient();
        FTPReply reply = ftp.connect(host);
        assertTrue(reply.isSuccess());

        // TODO -- handle failures
        reply = ftp.login("anonymous", "igv@broadinstitute.org");
        assertTrue(reply.isSuccess());

        reply = ftp.binary();
        assertTrue(reply.isSuccess());

        return ftp;

    }
}

