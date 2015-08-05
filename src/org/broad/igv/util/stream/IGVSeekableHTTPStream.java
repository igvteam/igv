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
import org.apache.log4j.Logger;
import org.broad.igv.ui.util.Packable;
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

        int attempts = 0;
        while(attempts < 3) {
            try {
                return _read(buffer, offset, len);
            } catch (java.net.SocketException e) {
                if(attempts < 3) {
                    attempts++;
                    log.error("Socket exception. Trying again.", e);
                }
                else {
                    throw e;
                }
            }
        }

        throw new RuntimeException("Reading " + url + " failed with unknown error.");  // Should be impossible to get here
    }

    public int _read(byte[] buffer, int offset, int len) throws IOException {

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
