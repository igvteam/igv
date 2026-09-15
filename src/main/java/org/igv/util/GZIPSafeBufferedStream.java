package org.igv.util;

import java.io.BufferedInputStream;
import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Objects;

/**
 * A BufferedInputStream whose available() returns 0 only at the end of the stream, blocking if necessary until a
 * byte can be read.
 * <p>
 * This works around GZIPInputStream stopping early on concatenated gzip members, e.g. bgzipped files, read from a
 * network stream.  After each member's trailer GZIPInputStream calls available() on its input to decide whether
 * another member follows, but network streams can report 0 before the end of the stream.  The available() check was
 * removed in JDK 23 (backported to 21.0.9) but restored in JDK 27 (backported to 25.0.5 and 26.0.2), where setting
 * the jdk.util.gzip.tryReadAheadAfterTrailer system property skips it.
 * <p>
 * See https://github.com/igvteam/igv/issues/1693 and https://github.com/samtools/htsjdk/issues/1691
 */
public class GZIPSafeBufferedStream extends BufferedInputStream {

    private final LookaheadInputStream lookahead;

    public GZIPSafeBufferedStream(InputStream in) {
        this(new LookaheadInputStream(in));
    }

    private GZIPSafeBufferedStream(LookaheadInputStream lookahead) {
        super(lookahead);
        this.lookahead = lookahead;
    }

    @Override
    public synchronized int available() throws IOException {
        int n = super.available();
        return n == 0 && lookahead.fetch() ? super.available() : n;
    }

    /**
     * Holds at most one byte read ahead from the wrapped stream, leaving the enclosing BufferedInputStream's buffer
     * and mark undisturbed.  Its own available() never blocks, as BufferedInputStream calls it while reading.
     */
    private static class LookaheadInputStream extends FilterInputStream {

        private static final int NONE = -2;

        private int next = NONE;   // NONE, -1 for end of stream, or the byte read ahead

        LookaheadInputStream(InputStream in) {
            super(in);
        }

        /**
         * Block until a byte has been read ahead or the end of the stream is reached.
         *
         * @return true if a byte is available
         */
        boolean fetch() throws IOException {
            if (next == NONE) next = in.read();
            return next >= 0;
        }

        @Override
        public int read() throws IOException {
            if (next == NONE) return in.read();
            int b = next;
            if (b >= 0) next = NONE;
            return b;
        }

        @Override
        public int read(byte[] b, int off, int len) throws IOException {
            Objects.checkFromIndexSize(off, len, b.length);
            if (next == NONE) return in.read(b, off, len);
            if (len == 0) return 0;
            if (next < 0) return -1;
            b[off] = (byte) next;
            next = NONE;
            // A byte has been read, don't block waiting for more
            int n = len > 1 && in.available() > 0 ? in.read(b, off + 1, len - 1) : 0;
            return 1 + Math.max(n, 0);
        }

        @Override
        public long skip(long n) throws IOException {
            if (next == NONE) return in.skip(n);
            if (n <= 0 || next < 0) return 0;
            next = NONE;
            return 1 + in.skip(n - 1);
        }

        @Override
        public int available() throws IOException {
            if (next == NONE) return in.available();
            if (next < 0) return 0;
            int n = in.available();
            return n == Integer.MAX_VALUE ? n : n + 1;
        }
    }
}
