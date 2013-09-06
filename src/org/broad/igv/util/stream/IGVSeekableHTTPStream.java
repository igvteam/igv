/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.util.stream;

import net.sf.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.tribble.util.URLHelper;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * TODO Get rid of this class
 * Temporary workaround to let us use IGVUrlHelper, think
 * it's what's causing the issue with reading TDF from web start over http
 */
public class IGVSeekableHTTPStream extends SeekableStream {

    static Logger log = Logger.getLogger(IGVSeekableHTTPStream.class);

    private long position = 0;
    private long contentLength = -1;

    private URLHelper helper;

    public IGVSeekableHTTPStream(final URL url) {

        this.helper = new IGVUrlHelper(url);
        log.debug("Getting content length for " + url);
        try {
            this.contentLength = this.helper.getContentLength();
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException(e.getMessage(), e);
        }

    }

    public long position() {
        return position;
    }

    public long length() {
        return contentLength;
    }

    @Override
    public long skip(long n) throws IOException {
        long bytesToSkip = Math.min(n, contentLength - position);
        position += bytesToSkip;
        return bytesToSkip;
    }

    public boolean eof() throws IOException {
        return contentLength > 0 && position >= contentLength;
    }

    public void seek(final long position) {
        this.position = position;
    }

    public int read(byte[] buffer, int offset, int len) throws IOException {

        String stats = "Offset="+offset+",len="+len+",buflen="+buffer.length;
        log.debug("Reading from " + getSource());
        log.debug("Stats: " + stats);
        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
            throw new IndexOutOfBoundsException(stats);
        }
        if (len == 0) {
            return 0;
        }

        InputStream is = null;
        int n = 0;
        try {

            long endRange = position + len - 1;
            // IF we know the total content length, limit the end range to that.
            if (contentLength > 0) {
                endRange = Math.min(endRange, contentLength);
            }
            log.debug("Trying to read range " + position + " to " + endRange);
            is = this.helper.openInputStreamForRange(position, endRange);

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
            log.debug("Read " + n  + " bytes, current position now " + position);
            return n;

        }

        catch (IOException e) {
            // THis is a bit of a hack, but its not clear how else to handle this.  If a byte range is specified
            // that goes past the end of the file the response code will be 416.  The MAC os translates this to
            // an IOException with the 416 code in the message.  Windows translates the error to an EOFException.
            //
            //  The BAM file iterator  uses the return value to detect end of file (specifically looks for n == 0).
            log.error(e.getMessage(), e);
            if (e.getMessage().contains("416") || (e instanceof EOFException)) {
                if (n == 0) {
                    return -1;
                } else {
                    position += n;
                    // As we are at EOF, the contentLength and position are by definition =
                    contentLength = position;
                    return n;
                }
            } else {
                throw e;
            }

        }

        finally {
            if (is != null) {
                is.close();
            }
        }
    }


    public void close() throws IOException {
        // Nothing to do
    }


    public int read() throws IOException {
        byte []tmp=new byte[1];
        read(tmp,0,1);
        return (int) tmp[0] & 0xFF;
    }

    @Override
    public String getSource() {
        return this.helper.getUrl().toExternalForm();
    }
}
