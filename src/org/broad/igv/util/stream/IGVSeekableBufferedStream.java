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

import com.google.common.primitives.Ints;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;

import java.io.IOException;

import static java.lang.System.arraycopy;

/**
 * A wrapper class to provide buffered read access to a SeekableStream.  Just wrapping such a stream with
 * a BufferedInputStream will not work as it does not support seeking.  In this implementation,
 * we attempt to reuse the buffer if there is overlap between the newly requested range and where
 * the buffer contains data for.
 */
public class IGVSeekableBufferedStream extends SeekableStream {

    private static Logger log = Logger.getLogger(IGVSeekableBufferedStream.class);

    public static final int DEFAULT_BUFFER_SIZE = 512000;

    final private int maxBufferSize;
    final SeekableStream wrappedStream;
    long position;
    /**
     * Can't make the assumption that length is a valid value.
     * May be -1, which means we don't know.
     */

    int markpos;
    int marklimit;

    byte[] buffer;
    long bufferStartPosition; // Position in file corresponding to start of buffer
    int bufferSize;

    public IGVSeekableBufferedStream(final SeekableStream stream, final int bsize) {
        this.maxBufferSize = bsize;
        this.wrappedStream = stream;
        this.position = 0;
        this.buffer = new byte[maxBufferSize];
        this.bufferStartPosition = -1;
        this.bufferSize = 0;

    }

    public IGVSeekableBufferedStream(final SeekableStream stream) {
        this(stream, DEFAULT_BUFFER_SIZE);
    }

    public long length() {
        return wrappedStream.length();
    }

    @Override
    public long skip(final long skipLength) throws IOException {
        long maxSkip = Long.MAX_VALUE;
        long length = wrappedStream.length();
        if (length >= 0) maxSkip = length - position - 1;
        long actualSkip = Math.min(maxSkip, skipLength);
        position += actualSkip;
        return actualSkip;
    }

    @Override
    public synchronized void reset() throws IOException {

        if (markpos < 0) {
            throw new IOException("Resetting to invalid mark");
        }
        position = markpos;

    }

    @Override
    public synchronized void mark(int readlimit) {
        this.markpos = (int) position;
        this.marklimit = readlimit;
    }

    @Override
    public boolean markSupported() {
        return true;
    }

    public void seek(final long position) throws IOException {
        this.position = position;
    }

    public void close() throws IOException {
        wrappedStream.close();
    }

    public boolean eof() throws IOException {
        long length = wrappedStream.length();
        return length >= 0 && position >= length;
    }

    @Override
    public String getSource() {
        return wrappedStream.getSource();
    }

    @Override
    public long position() throws IOException {
        return position;
    }

    /**
     * Return true iff the buffer needs to be refilled for the given
     * amount of data requested
     *
     * @param len Number of bytes from {@code position} one plans on reading
     * @return
     */
    private boolean needFillBuffer(int len) {
        return bufferSize == 0 || position < bufferStartPosition || (position + len) > bufferStartPosition + bufferSize;
    }


    public int read() throws IOException {

        if (needFillBuffer(1)) {
            fillBuffer();
        }

        int offset = (int) (position - bufferStartPosition);
        int b = buffer[offset];
        position++;
        return b;
    }

    public int read(final byte[] b, final int off, final int len) throws IOException {

        long length = wrappedStream.length();
        if (length >= 0 && position >= length) return -1;

        if (len > maxBufferSize) {
            // Buffering not useful here.  Don't bother trying to use any (possible) overlapping buffer contents
            wrappedStream.seek(position);
            int nBytes = wrappedStream.read(b, off, len);
            position += nBytes;
            return nBytes;
        } else {
            // Requested range is not contained within buffer.
            if (needFillBuffer(len)) {
                fillBuffer();
            }

            int bufferOffset = (int) (position - bufferStartPosition);
            int bytesCopied = Math.min(len, bufferSize - bufferOffset);
            arraycopy(buffer, bufferOffset, b, off, bytesCopied);
            position += bytesCopied;
            return bytesCopied;
        }
    }

    private void fillBuffer() throws IOException {

        long longRem = maxBufferSize;
        long length = wrappedStream.length();

        if (length >= 0) longRem = Math.min((long) maxBufferSize, length - position);

        //This shouldn't actually be necessary as long as maxBufferSize is
        //an int, but we leave it here to stress the fact that
        //we need to watch for overflow
        int bytesRemaining = Ints.saturatedCast(longRem);
        //Number of bytes to skip at beginning when reading from stream later
        int toSkip = 0;
        //Number of bytes known to be stored in the buffer, which are valid
        int tmpBufferSize = 0;

        long bufferEnd = bufferStartPosition + bufferSize;
        long requiredEnd = position + bytesRemaining;

        if (bufferStartPosition > position && bufferEnd < requiredEnd) {
            // Buffer is in the middle of the required range, don't bother trying to reuse
        } else {
            if (position < bufferEnd && requiredEnd > bufferEnd) {
                // Beginning of buffer data is useless, want to save
                // some though
                // Fill request:           xxxxxxxx...
                // Buffer data:    xxxxxxxxxxxxx
                int szOverlap = (int) (bufferEnd - position);
                try {
                    arraycopy(buffer, bufferSize - szOverlap, buffer, 0, szOverlap);
                } catch (Exception e) {
                    e.printStackTrace();
                }

                //Skip the first bytes that we already had
                toSkip = szOverlap;
                tmpBufferSize += szOverlap;
                bytesRemaining -= szOverlap;
            } else if (position < bufferStartPosition && requiredEnd > bufferStartPosition && requiredEnd < bufferEnd) {
                //Gap between position and buffer start, but some overlap
                //We require that the buffer contain data all the way to the end,
                //because dealing with the case of writing buffered data to the middle
                //is too complicated and not likely to occur in practice. When it does,
                //we just re-read everything.
                // Fill request: xxxxxxxxxxxx...
                // Buffer data:      xxxxxxxxxxx...
                int szOverlap = (int) (requiredEnd - bufferStartPosition);
                int destPos = (int) (bufferStartPosition - position);
                arraycopy(buffer, 0, buffer, destPos, szOverlap);

                //Don't skip any bytes, just trim from the number requested
                bytesRemaining -= szOverlap;
                tmpBufferSize += szOverlap;
            }
        }

        if (bytesRemaining > 0) {
            int curOffset = toSkip;
            wrappedStream.seek(position + toSkip);
            while (bytesRemaining > 0) {
                int count = wrappedStream.read(buffer, curOffset, bytesRemaining);
                if (count < 0) {
                    break;  // EOF.
                }
                curOffset += count;
                bytesRemaining -= count;
                tmpBufferSize += count;
            }
            bufferStartPosition = position;
            bufferSize = tmpBufferSize;
        }
    }
}