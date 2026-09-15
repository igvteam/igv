package org.igv.util;

import org.junit.Test;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Random;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import static org.junit.Assert.*;

public class GZIPSafeBufferedStreamTest {

    /**
     * Simulates a network stream: returns at most one chunk per read and always reports 0 bytes available.
     */
    static class ChunkedStream extends InputStream {

        private final byte[] data;
        private final TreeSet<Integer> chunkEnds = new TreeSet<>();
        private int pos = 0;

        ChunkedStream(byte[] data, Iterable<Integer> chunkEnds) {
            this.data = data;
            chunkEnds.forEach(this.chunkEnds::add);
        }

        ChunkedStream(byte[] data, int chunkSize) {
            this.data = data;
            for (int i = chunkSize; i < data.length; i += chunkSize) chunkEnds.add(i);
        }

        @Override
        public int read() {
            return pos < data.length ? data[pos++] & 0xff : -1;
        }

        @Override
        public int read(byte[] b, int off, int len) {
            if (len == 0) return 0;
            if (pos >= data.length) return -1;
            Integer end = chunkEnds.higher(pos);
            int n = Math.min(len, (end == null ? data.length : end) - pos);
            System.arraycopy(data, pos, b, off, n);
            pos += n;
            return n;
        }

        @Override
        public int available() {
            return 0;
        }
    }

    @Test
    public void testAvailableUntilEndOfStream() throws IOException {

        byte[] data = new byte[1000];
        new Random(1).nextBytes(data);

        GZIPSafeBufferedStream stream = new GZIPSafeBufferedStream(new ChunkedStream(data, 7));
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        byte[] buf = new byte[5];
        while (out.size() < data.length) {
            assertTrue(stream.available() > 0);
            if (out.size() % 2 == 0) {
                out.write(stream.read());
            } else {
                int n = stream.read(buf);
                out.write(buf, 0, n);
            }
        }
        assertEquals(0, stream.available());
        assertEquals(-1, stream.read());
        assertArrayEquals(data, out.toByteArray());
    }

    /**
     * ParsingUtils.isGzip marks the stream, reads the 2 byte signature, and resets.  A first chunk of exactly 2 bytes
     * must not lose the signature.
     */
    @Test
    public void testMarkResetPreserved() throws IOException {

        String text = "chr1\t100\t200\tfeature\n".repeat(100);
        byte[] gzipped = gzip(text);

        GZIPSafeBufferedStream stream = new GZIPSafeBufferedStream(new ChunkedStream(gzipped, java.util.List.of(2)));
        stream.mark(2);
        assertEquals(2, stream.read(new byte[2]));
        assertTrue(stream.available() > 0);
        stream.reset();

        assertTrue(ParsingUtils.isGzip(stream));
        stream.reset();
        assertEquals(text, new String(new GZIPInputStream(stream).readAllBytes(), StandardCharsets.UTF_8));
    }

    @Test
    public void testSkipAfterLookahead() throws IOException {

        byte[] data = new byte[100];
        for (int i = 0; i < data.length; i++) data[i] = (byte) i;

        GZIPSafeBufferedStream stream = new GZIPSafeBufferedStream(new ChunkedStream(data, 10));
        assertEquals(10, stream.read(new byte[10]));
        assertTrue(stream.available() > 0);
        assertEquals(5, stream.skip(5));
        assertEquals(15, stream.read());
    }

    /**
     * Concatenated gzip members with chunks ending at member boundaries and within member trailers, read directly
     * (as ParsingUtils.openInputStreamGZ does) and through an intermediate buffer (as htsjdk does).
     */
    @Test
    public void testConcatenatedMembers() throws IOException {

        ByteArrayOutputStream concatenated = new ByteArrayOutputStream();
        StringBuilder expected = new StringBuilder();
        TreeSet<Integer> chunkEnds = new TreeSet<>();
        for (int m = 0; m < 5; m++) {
            String text = ("chr1\t" + m + "\t200\tfeature\n").repeat(1000 + m * 37);
            expected.append(text);
            concatenated.write(gzip(text));
            chunkEnds.add(concatenated.size() - 3);
            chunkEnds.add(concatenated.size());
        }
        byte[] data = concatenated.toByteArray();

        InputStream direct = new GZIPInputStream(new GZIPSafeBufferedStream(new ChunkedStream(data, chunkEnds)));
        assertEquals(expected.toString(), new String(direct.readAllBytes(), StandardCharsets.UTF_8));

        InputStream buffered = new GZIPInputStream(new BufferedInputStream(new GZIPSafeBufferedStream(new ChunkedStream(data, chunkEnds)), 512000));
        assertEquals(expected.toString(), new String(buffered.readAllBytes(), StandardCharsets.UTF_8));
    }

    private static byte[] gzip(String text) throws IOException {
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        try (OutputStream out = new GZIPOutputStream(bos)) {
            out.write(text.getBytes(StandardCharsets.UTF_8));
        }
        return bos.toByteArray();
    }
}
