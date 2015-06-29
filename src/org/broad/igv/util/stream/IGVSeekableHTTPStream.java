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

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.util.HttpUtils;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 */
public class IGVSeekableHTTPStream extends SeekableStream {

    static Logger log = Logger.getLogger(IGVSeekableHTTPStream.class);

    private long position = 0;
    private URL url;
    long contentLength = -1;                      // Not set

    public IGVSeekableHTTPStream(final URL url) {
        this.url = url;
    }

    public long position() {
        return position;
    }

    public long length() {
        return contentLength;
    }

    @Override
    public long skip(long n) throws IOException {
        long bytesToSkip = contentLength < 0 ? n : Math.min(n, contentLength - position);
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

        String stats = "Offset=" + offset + ",len=" + len + ",buflen=" + buffer.length;
        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
            throw new IndexOutOfBoundsException(stats);
        }
        if (len == 0) {
            return 0;
        }

        InputStream is = null;
        int n = 0;
        try {

            if (contentLength > 0 && position >= contentLength) {
                return -1;  // EOF
            }

            long endRange = position + len - 1;
            // IF we know the total content length, limit the end range to that.
            if (contentLength > 0) {
                endRange = Math.min(endRange, contentLength);
            }
            if (log.isTraceEnabled()) {
                log.trace("Trying to read range " + position + " to " + endRange);
            }
            is = openInputStreamForRange(position, endRange);

            while (n < len) {
                int count = robustRead(buffer, offset + n, len - n, is, 0);
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

        } catch (IOException e) {
            // THis is a bit of a hack, but its not clear how else to handle this.  If a byte range is specified
            // that goes past the end of the file the response code will be 416.  The MAC os translates this to
            // an IOException with the 416 code in the message.  Windows translates the error to an EOFException.
            //
            //  The BAM file iterator  uses the return value to detect end of file (specifically looks for n == 0).
            if (e.getMessage().contains("416") || (e instanceof EOFException)) {
                if (n == 0) {
                    contentLength = position;
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

        } finally {
            if (is != null) {
                is.close();
            }
        }
    }

    public int robustRead(byte[] buffer, int offset, int len, InputStream is, int attempts) throws IOException {
        try {
            return is.read(buffer, offset, len);
        } catch (java.net.SocketException e) {
            if (attempts < 4) {
                return robustRead(buffer, offset, len, is, attempts + 1);
            } else {
                throw e;
            }
        }
    }


    public void close() throws IOException {
        // Nothing to do
    }


    public int read() throws IOException {
        byte[] tmp = new byte[1];
        read(tmp, 0, 1);
        return (int) tmp[0] & 0xFF;
    }

    public InputStream openInputStreamForRange(long start, long end) throws IOException {

        String byteRange = "bytes=" + start + "-" + end;
        Map<String, String> params = new HashMap();
        params.put("Range", byteRange);
        //URL url = addStartEndQueryString(this.url, start, end);

        HttpURLConnection conn = HttpUtils.getInstance().openConnection(url, params);
        try {
            int contentLength = conn.getContentLength();
            if (contentLength > 0 && (contentLength != (end - start + 1))) {
                // We're at EOF, record
                this.contentLength = this.position + contentLength;
            }
        } catch (Exception e) {
            log.error("Error determining content length", e);
        }

        try {
            InputStream input = conn.getInputStream();
            return input;
        } catch (IOException e) {
            HttpUtils.getInstance().readErrorStream(conn);  // Consume content
            throw e;
        }
    }


    /**
     * Add query parameters which should more properly be in Range header field
     * to query string
     *
     * @param start start byte
     * @param end   end byte
     * @throws java.net.MalformedURLException
     */
    static URL addStartEndQueryString(URL url, long start, long end) throws MalformedURLException {

        String queryString = url.getQuery();
        if (queryString == null) {
            return new URL(url.toExternalForm() + "?start=" + start + "&end=" + end);
        } else {
            String newQueryString = queryString + "&start=" + start + "&end" + end;
            return new URL(url.toExternalForm().replace(queryString, newQueryString));
        }
    }

    @Override
    public String getSource() {
        return url.toExternalForm();
    }


    public static void main(String[] args) throws IOException {

        IGVSeekableHTTPStream stream = new IGVSeekableHTTPStream(new URL("http://localhost/igv-web/test/data/misc/BufferedReaderTest.bin"));

        byte[] buffer = new byte[1000];

        stream.read(buffer, 0, 1000);

        System.out.println("Done");

    }
}
