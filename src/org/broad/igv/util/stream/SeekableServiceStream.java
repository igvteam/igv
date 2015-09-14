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
import org.broad.igv.util.HttpUtils;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

/**
 * A SeekableStream implementation for the "range" webservice.  The purpose of this class is to serve range-byte
 * requests to clients who are unable to use the http header for this.
 * <p/>
 * /xchip/igv/data/public/annotations/seq/hg18/chr1.txt
 */
public class SeekableServiceStream extends SeekableStream {

    static Logger log = Logger.getLogger(SeekableServiceStream.class);

    private static final String WEBSERVICE_URL = "http://www.broadinstitute.org/webservices/igv/range";



    private long position = 0;
    private long contentLength = Long.MAX_VALUE;
    private URL  wrappedURL;

    public SeekableServiceStream(URL url) {
        this.wrappedURL = url;
    }

    public long length() {
        if(contentLength == Long.MAX_VALUE) {
            try {
                contentLength = HttpUtils.getInstance().getContentLength(wrappedURL);
            } catch (IOException e) {
                log.error("Error fetching content length for: " + wrappedURL, e);
            }
        }
        return contentLength;
    }

    public boolean eof() throws IOException {
        return position >= contentLength;
    }

    @Override
    public String getSource() {
        return this.wrappedURL.toExternalForm();
    }

    public void seek(long position) {
        this.position = position;
    }

    public long position() {
        return position;
    }

    @Override
    public long skip(long n) throws IOException {
        position += n;
        return n;
    }

    public int read(byte[] buffer, int offset, int length) throws IOException {
        if (offset < 0 || length < 0 || (offset + length) > buffer.length) {
            throw new IndexOutOfBoundsException();
        }

        InputStream is = null;

        URL url = new URL(WEBSERVICE_URL + "?file=" + wrappedURL.toExternalForm() + "&position=" + position + "&length=" + length);

        int n = 0;
        try {

            is = HttpUtils.getInstance().openConnectionStream(url);

            while (n < length) {
                int count = is.read(buffer, offset + n, length - n);
                if (count < 0) {
                    return n;
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
                return n;
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


    public byte[] readBytes(long position, int nBytes) throws IOException {
        this.position = position;
        byte[] buffer = new byte[nBytes];
        read(buffer, 0, nBytes);
        return buffer;
    }

    public int read() throws IOException {
        throw new UnsupportedOperationException("read() is not supported on SeekableServiceStream.  Must read in blocks.");
    }
}