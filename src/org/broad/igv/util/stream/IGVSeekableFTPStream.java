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

package org.broad.igv.util.stream;

import org.apache.log4j.Logger;
import org.broad.igv.util.UserPasswordInputImpl;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.ftp.FTPClient;
import org.broad.tribble.util.ftp.FTPReply;
import org.broad.tribble.util.ftp.FTPUtils;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import static junit.framework.Assert.assertTrue;

/**
 * @author jrobinso
 * @date Oct 27, 2010
 */
public class IGVSeekableFTPStream extends SeekableStream {

    private static Logger log = Logger.getLogger(IGVSeekableFTPStream.class);

    private long position = 0;
    private String host;
    private String path;
    private String userInfo;
    FTPClient ftp = null;

    public IGVSeekableFTPStream(URL url) throws IOException {
        this.userInfo = url.getUserInfo();
        this.host = url.getHost();
        this.path = url.getPath();
        ftp = FTPUtils.connect(host, userInfo, new UserPasswordInputImpl());
    }

    public void seek(long position) {
        this.position = position;
    }

    public long position() {
        return position;
    }

    public boolean eof() throws IOException {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public long length() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }


    public long skip(long n) throws IOException {
        long bytesToSkip = n;
        position += bytesToSkip;
        ftp.setRestPosition(position);
        return bytesToSkip;
    }

    public int read(byte[] buffer, int offset, int len) throws IOException {

        if (ftp == null) {
            ftp = FTPUtils.connect(host, userInfo, new UserPasswordInputImpl());
        }

        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
            throw new IndexOutOfBoundsException();
        }

        if (len == 0) {
            return 0;
        }

        int n = 0;
        try {

            FTPReply reply = ftp.pasv();
            assertTrue(reply.isSuccess());

            if (position > 0) ftp.setRestPosition(position);

            reply = ftp.retr(path);
            //assertTrue(reply.getCode() == 150);

            InputStream is = ftp.getDataStream();

            while (n < len) {
                int count = is.read(buffer, offset + n, len - n);
                if (count < 0) {
                    if (n == 0) {
                        return -1;
                    } else {
                        break;
                    }
                }
                n += count;
            }

            position += n;

            return n;
        }

        catch (EOFException e) {
            if (n < 0) {
                return -1;
            } else {
                position += n;
                return n;
            }
        }

        finally {
            // ALWAYS close ftp connection,  this is more robust than trying to resue them,
            // and we don't want open connections hanging about
            ftp.disconnect();
            ftp = null;
        }
    }


    private void reconnect() throws IOException {
        if (ftp != null) {
            ftp.disconnect();
        }
        ftp = FTPUtils.connect(host, userInfo, new UserPasswordInputImpl());
    }


    public void close() throws IOException {
        log.info("close");
        if (ftp != null) {
            ftp.disconnect();
            ftp = null;
        }
    }


    public int read() throws IOException {
        throw new UnsupportedOperationException("read() is not supported on SeekableHTTPStream.  Must read in blocks.");
    }

}
