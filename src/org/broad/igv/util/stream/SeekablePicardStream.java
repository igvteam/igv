package org.broad.igv.util.stream;


import org.broad.tribble.util.SeekableHTTPStream;

import java.io.IOException;
import java.net.URL;

/**
 * Unfortunately the seekable stream classes exist for both Tribble and Picard, and we need both.  This class
 * is for use with Picard (SAM) and delegates all the work to a helper.
 *
 * @author jrobinso
 * @date Dec 23, 2011
 */

public class SeekablePicardStream extends net.sf.samtools.util.SeekableStream  {

    org.broad.tribble.util.SeekableStream tribbleStream;
    String source;

    public SeekablePicardStream(org.broad.tribble.util.SeekableStream tribbleStream, String source) {
        this.tribbleStream = tribbleStream;
        this.source = source;
    }

    @Override
    public long length() {
       return tribbleStream.length();
    }

    @Override
    public void seek(long l) throws IOException {
        tribbleStream.seek(l);
    }

    @Override
    public int read() throws IOException {
        return tribbleStream.read();
    }

    @Override
    public int read(byte[] bytes, int offset, int length) throws IOException {
        return tribbleStream.read(bytes,  offset,  length);
    }

    @Override
    public void close() throws IOException {
        tribbleStream.close();
    }

    @Override
    public boolean eof() throws IOException {
        return tribbleStream.eof();
    }

    @Override
    public String getSource() {
        return source;
    }
}
