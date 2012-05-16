package org.broad.igv.util.stream;


import net.sf.samtools.util.SeekableStream;

import java.io.BufferedInputStream;
import java.io.IOException;

/**
 * Unfortunately the seekable stream classes exist for both Tribble and Picard, and we need both.  This class
 * is for use with Picard (SAM) and delegates all the work to a helper.
 *
 * @author jrobinso
 * @date Dec 23, 2011
 */

public class SeekablePicardStream extends net.sf.samtools.util.SeekableStream {

    public static final int DEFAULT_BUFFER_SIZE = 1024000;

    final private int bufferSize;
    org.broad.tribble.util.SeekableStream wrappedStream;
    String source;
    BufferedInputStream bufferedStream;
    long position;


    public SeekablePicardStream(org.broad.tribble.util.SeekableStream tribbleStream, String source) {

        this.wrappedStream = tribbleStream;
        this.source = source;
        this.position = 0;
        this.bufferSize = DEFAULT_BUFFER_SIZE;
        bufferedStream = new BufferedInputStream(wrappedStream, bufferSize);
    }

    @Override
    public long length() {
        return wrappedStream.length();
    }

    @Override
    public void seek(long position) throws IOException {
        this.position = position;
        wrappedStream.seek(position);
        bufferedStream = new BufferedInputStream(wrappedStream, bufferSize);
    }

    @Override
    public int read() throws IOException {
        int b = bufferedStream.read();
        position++;
        return b;
    }


    @Override
    public int read(byte[] buffer, int offset, int length) throws IOException {
        int nBytesRead = bufferedStream.read(buffer, offset, length);
        if (nBytesRead > 0) {
            position += nBytesRead;
        }
        return nBytesRead;
    }


    @Override
    public void close() throws IOException {
        wrappedStream.close();
    }

    @Override
    public boolean eof() throws IOException {
        return wrappedStream.eof();
    }

    @Override
    public String getSource() {
        return source;
    }
}
