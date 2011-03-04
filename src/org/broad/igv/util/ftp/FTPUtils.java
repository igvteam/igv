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


import org.apache.log4j.Logger;
import org.broad.igv.ui.UserPasswordPane;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.Hashtable;
import java.util.Map;


/**
 * @author jrobinso
 * @date Aug 31, 2010
 */
public class FTPUtils {

    private static Logger log = Logger.getLogger(FTPUtils.class);

    static Map<String, String> userCredentials = new Hashtable();

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
     * ftp://virothao:satio-ngeca-62@shiraz.well.ox.ac.uk/110210_SN228_0119_A702MDABXX/WTCHG_11474/WTCHG_11474_101.bam
     * ftp://shiraz.well.ox.ac.uk/110210_SN228_0119_A702MDABXX/WTCHG_11474/WTCHG_11474_101.bam
     * virothao:satio-ngeca-62@
     * <p/>
     * Connect to an FTP server anonymously
     *
     * @param host
     * @param userInfo
     * @return
     * @throws IOException
     */
    public static FTPClient connect(String host, String userInfo) throws IOException {

        log.info("connect");

        FTPClient ftp = new FTPClient();
        FTPReply reply = ftp.connect(host);
        if (!reply.isSuccess()) {
            throw new RuntimeException("Could not connect to " + host);
        }

        String user = "anonymous";
        String password = "igv@broadinstitute.org";

        if (userInfo == null) {
            userInfo = userCredentials.get(host);
        }
        if (userInfo != null) {
            String[] tmp = userInfo.split(":");
            user = tmp[0];
            if (tmp.length > 1) {
                password = tmp[1];
            }
        }

        reply = ftp.login(user, password);
        if (!reply.isSuccess()) {

            UserPasswordPane pane = new UserPasswordPane(host);
            boolean success = false;
            while (!success) {
                if (pane.showDialog()) {
                    user = pane.getUser();
                    password = pane.getPassword();
                    reply = ftp.login(user, password);
                    success = reply.isSuccess();
                } else {
                    // canceled
                    break;
                }

            }
            if (success) {
                userInfo = user + ":" + password;
                userCredentials.put(host, userInfo);
            } else {
                throw new RuntimeException("Login failure for host: " + host);
            }
        }

        reply = ftp.binary();
        if (!(reply.isSuccess())) {
            throw new RuntimeException("Could not set binary mode on host: " + host);
        }

        log.info("Connected to " + host);

        return ftp;

    }

}

