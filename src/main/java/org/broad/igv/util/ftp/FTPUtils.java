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


import htsjdk.samtools.SAMException;
import htsjdk.samtools.seekablestream.UserPasswordInput;
import htsjdk.samtools.util.ftp.FTPClient;
import htsjdk.samtools.util.ftp.FTPReply;


import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.HashMap;
import java.util.Map;


/**
 * @author jrobinso
 * @date Aug 31, 2010
 */
public class FTPUtils {

    static Map<String, String> userCredentials = new HashMap<String, String>();

    static int TIMEOUT = 10000;

    public static boolean resourceAvailable(URL url) {
        InputStream is = null;
        try {

            return getContentLength(url) > 0;

        } catch (IOException e) {
            return false;
        }
        // NOTE-- DO NOT TRY TO CLOSE STREAM.  IT WILL HANG.
    }

    public static long getContentLength(URL url) throws IOException {
        FTPClient ftp = null;
        try {
            ftp = FTPUtils.connect(url.getHost(), url.getUserInfo(), null);
            String sizeString = ftp.executeCommand("size " + url.getPath()).getReplyString();
            return Integer.parseInt(sizeString);
        } catch (Exception e) {
            return -1 ;
        }
        finally {
            if(ftp != null) {
                ftp.disconnect();
            }
        }

    }


    /**
     * Connect to an FTP server
     *
     * @param host
     * @param userInfo
     * @param userPasswordInput Dialog with which a user can enter credentials, if login fails
     * @return
     * @throws IOException
     */
    public static synchronized FTPClient connect(String host, String userInfo, UserPasswordInput userPasswordInput) throws IOException {

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
            if (userPasswordInput == null) {
                throw new RuntimeException("Login failure for host: " + host);
            } else {
                userPasswordInput.setHost(host);
                boolean success = false;
                while (!success) {
                    if (userPasswordInput.showDialog()) {
                        user = userPasswordInput.getUser();
                        password = userPasswordInput.getPassword();
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
        }

        reply = ftp.binary();
        if (!(reply.isSuccess())) {
            throw new RuntimeException("Could not set binary mode on host: " + host);
        }

        return ftp;

    }

}

