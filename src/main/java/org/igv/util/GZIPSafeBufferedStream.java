package org.igv.util;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Work around for a bug in GZIPInputStream (fixed in JDK 23) when reading concatenated gzip members, e.g. bgzipped
 * files, from a network stream.  GZIPInputStream uses available() to decide if another member follows, but network
 * streams can report 0 bytes available before the end of the stream.  This class keeps available() > 0 until the
 * end of the stream by pre-filling the buffer.
 * <p>
 * See https://github.com/igvteam/igv/issues/1693 and https://github.com/samtools/htsjdk/issues/1691
 */
public class GZIPSafeBufferedStream extends BufferedInputStream {

    public GZIPSafeBufferedStream(InputStream in) {
        super(in);
    }


    @Override
    public int read(byte[] b, int off, int len) throws IOException {
        int ret = super.read(b, off, len);
        if (available() == 0) {
            mark(26);
            read();
            reset();
        }
        return ret;
    }


}
