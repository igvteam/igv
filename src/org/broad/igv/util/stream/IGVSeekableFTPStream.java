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

package org.broad.igv.util.stream;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.ftp.FTPClient;
import htsjdk.samtools.util.ftp.FTPReply;
import org.broad.igv.util.ftp.FTPUtils;
import org.apache.log4j.Logger;
import org.broad.igv.util.UserPasswordInputImpl;

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
    private String source;

    private long length = -1;

    public IGVSeekableFTPStream(URL url) throws IOException {
        this.source = url.toExternalForm();
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
        return this.position() >= this.length();
    }

    @Override
    public String getSource() {
        return this.source;
    }

    public long length() {
        if(this.length < 0){
            this.length = getLength();
        }
        return this.length;
    }

    private long getLength(){
        try {
            FTPReply reply = ftp.size(path);
            if (reply.isSuccess()) {
                return Long.parseLong(reply.getReplyString());
            }
        } catch (IOException e) {
            log.error("Error getting length. " + e.getMessage(), e);
        }
        return -1;
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
        log.debug("close");
        if (ftp != null) {
            ftp.disconnect();
            ftp = null;
        }
    }


    public int read() throws IOException {
        throw new UnsupportedOperationException("read() is not supported on SeekableHTTPStream.  Must read in blocks.");
    }

}
