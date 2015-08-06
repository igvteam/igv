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
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableHTTPStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Assume;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;
import java.net.URL;
import java.util.Arrays;
import java.util.Random;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class IGVSeekableBufferedStreamTest extends AbstractHeadlessTest{

    //    private final File BAM_INDEX_FILE = new File("testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam.bai");
    private final File BAM_FILE = new File(TestUtils.DATA_DIR + "/samtools/index_test.bam");
    private final String BAM_URL_STRING = "http://www.broadinstitute.org/igvdata/test/data/bam/index_test.bam";
    private static File megabyteZerosFile = new File(TestUtils.DATA_DIR + "/samtools/megabyteZeros.dat");

    static int expectedFileSize = 20000;
    static byte[] expectedBytes;
    static File testFile;

    final int streamBufferSize = 100;
    final int halfStreamBufferSize = streamBufferSize / 2;
    IGVSeekableBufferedStream testBuffStream;


    @BeforeClass
    public static void setUpClass() throws Exception {
        createTestFile(expectedFileSize);
    }

    @Before
    public void setUp() throws Exception{
        testBuffStream = getTestStream(streamBufferSize);
    }


    /**
     * Test jumping around reading individual bytes
     * from a typical file stream
     *
     * @throws Exception
     */
    @Test
    public void testSeekReadStandard() throws Exception {
        tstSeekRead(expectedBytes, new SeekableFileStream(testFile), 50);
    }

    /**
     * Test jumping around reading individual bytes
     * from a stream of unknown length
     *
     * @throws Exception
     */
    @Test
    public void testSeekReadUnknownLength() throws Exception {
        tstSeekRead(expectedBytes, new UnknownLengthSeekableStream(testFile), 50);
    }

    /**
     * Test reading individual bytes
     *
     * @throws Exception
     */
    public void tstSeekRead(byte[] expectedBytes, SeekableStream wrappedStream, int bufferSize) throws Exception {

        final int fileSize = expectedBytes.length;

        IGVSeekableBufferedStream bufferedStream = new IGVSeekableBufferedStream(wrappedStream, bufferSize);

        // Somewhere in the middle
        int pos = 700;
        bufferedStream.seek(pos);
        byte b = (byte) bufferedStream.read();
        assertEquals(expectedBytes[pos], b);

        // near the end
        pos = fileSize - 100;
        bufferedStream.seek(pos);
        b = (byte) bufferedStream.read();
        assertEquals(expectedBytes[pos], b);

        // the beginning
        pos = 0;
        bufferedStream.seek(pos);
        b = (byte) bufferedStream.read();
        assertEquals(expectedBytes[pos], b);

        bufferedStream.close();
    }


    private IGVSeekableBufferedStream getTestStream(int streamBufferSize) throws FileNotFoundException{
       return new IGVSeekableBufferedStream(new SeekableFileStream(testFile), streamBufferSize);
    }
    /**
     * Test reading byte arrays
     *
     * @throws Exception
     */
    @Test
    public void testReadBuffer_noOverlap() throws Exception {

        byte[] buffer = new byte[halfStreamBufferSize];

        //Read at the beginning
        tstSeekReadBuffer(testBuffStream, expectedBytes, 0, buffer);

        //Jump to the end
        int pos = expectedFileSize - streamBufferSize;
        tstSeekReadBuffer(testBuffStream, expectedBytes, pos, buffer);

        // Somewhere in the middle
        pos = 700;
        tstSeekReadBuffer(testBuffStream, expectedBytes, pos, buffer);
    }


    @Test
    public void testReadBuffer_read_end() throws Exception{
        // Overlap with where the buffer has the beginning data, but not the end
        testBuffStream.seek(10000);
        testBuffStream.read();

        int pos = 10000 + halfStreamBufferSize;
        byte[] buffer = new byte[halfStreamBufferSize + 5];
        tstSeekReadBuffer(testBuffStream, expectedBytes, pos, buffer);
    }

    @Test
    public void testReadBuffer_missing_start_end() throws Exception{
        // Overlap with buffer missing data at the beginning and end
        // of the requested range.
        testBuffStream.seek(5000);
        testBuffStream.read();

        int pos = 5000 + halfStreamBufferSize;
        byte[] buffer = new byte[streamBufferSize];
        tstSeekReadBuffer(testBuffStream, expectedBytes, pos, buffer);
    }

    @Test
    public void testReadBuffer_read_start() throws Exception{
        // Overlap with stream missing data at the start of requested range
        testBuffStream.seek(5000);
        testBuffStream.read();

        int pos = 5000 - halfStreamBufferSize;
        byte[] buffer = new byte[halfStreamBufferSize - 5];
        tstSeekReadBuffer(testBuffStream, expectedBytes, pos, buffer);
    }

    @Test
    public void testReadBuffer_overlap_both() throws Exception{
        // Overlap both ends
        testBuffStream.seek(7000);
        testBuffStream.read();

        byte[] buffer = new byte[1000];
        int pos = 7000 - streamBufferSize - 25;
        tstSeekReadBuffer(testBuffStream, expectedBytes, pos, buffer);
    }

    @Test
    public void testReadBuffer_contained() throws Exception {
        // Completely contained
        testBuffStream.seek(3000);
        testBuffStream.read();

        int tstSize = streamBufferSize / 10;
        byte[] buffer = new byte[tstSize];
        int pos = 3000 + tstSize;
        tstSeekReadBuffer(testBuffStream, expectedBytes, pos, buffer);
    }


    private void tstSeekReadBuffer(IGVSeekableBufferedStream bufferedStream, byte[] allExpectedBytes, int pos, byte[] buffer) throws Exception{
        bufferedStream.seek(pos);
        bufferedStream.readFully(buffer);
        byte[] expected = Arrays.copyOfRange(allExpectedBytes, pos, pos + buffer.length);
        assertArraysEqual(expected, buffer);
    }

    private static void createTestFile(int length) throws Exception {

        expectedBytes = new byte[length];

        testFile = new File("SeekableBuffered_test.bin");
        testFile.deleteOnExit();

        int range = Byte.MAX_VALUE - Byte.MIN_VALUE;
        Random rand = new Random(19534);
        for (int i = 0; i < length; i++) {
            expectedBytes[i] = (byte) (Byte.MIN_VALUE + (int) (rand.nextDouble() * range));
        }

        OutputStream os = null;
        try {
            os = new BufferedOutputStream(new FileOutputStream(testFile));
            os.write(expectedBytes);
        } finally {
            if (os != null) os.close();
        }

    }

    /**
     * Test reading across a buffer boundary (buffer size is 512000).   The test first reads a range of
     * bytes using an unbuffered stream file stream,  then compares this to results from a buffered http stream.
     *
     * @throws IOException
     */
    @Test
    public void testRandomRead() throws IOException {

        int startPosition = 500000;
        int length = 50000;

        byte[] buffer1 = new byte[length];
        SeekableStream unBufferedStream = new SeekableFileStream(BAM_FILE);
        unBufferedStream.seek(startPosition);
        int bytesRead = unBufferedStream.read(buffer1, 0, length);
        assertEquals(length, bytesRead);

        byte[] buffer2 = new byte[length];
        SeekableStream bufferedStream = new IGVSeekableBufferedStream(new SeekableHTTPStream(new URL(BAM_URL_STRING)));
        bufferedStream.seek(startPosition);
        bytesRead = bufferedStream.read(buffer2, 0, length);
        assertEquals(length, bytesRead);

        assertArraysEqual(buffer1, buffer2);
    }


    /**
     * Test an attempt to read past the end of the file.  The test file is 594,149 bytes in length.  The test
     * attempts to read a 1000 byte block starting at position 594000.  A correct result would return 149 bytes.
     *
     * @throws IOException
     */
    @Test
    public void testEOF() throws IOException {

        int remainder = 149;
        long fileLength = BAM_FILE.length();
        long startPosition = fileLength - remainder;
        int length = 1000;


        byte[] buffer = new byte[length];
        SeekableStream bufferedStream = new IGVSeekableBufferedStream(new SeekableHTTPStream(new URL(BAM_URL_STRING)));
        bufferedStream.seek(startPosition);
        int bytesRead = bufferedStream.read(buffer, 0, length);
        assertEquals(remainder, bytesRead);

        // Subsequent reads should return -1
        bytesRead = bufferedStream.read(buffer, 0, length);
        assertEquals(-1, bytesRead);
    }

    @Test
    public void testSkip() throws IOException {
        final int[] BUFFER_SIZES = new int[]{8, 96, 1024, 8*1024, 16*1024, 96*1024, 48*1024};

        for (final int bufferSize : BUFFER_SIZES) {
            final IGVSeekableBufferedStream in1 = new IGVSeekableBufferedStream(new SeekableFileStream(BAM_FILE), bufferSize);
            final IGVSeekableBufferedStream in2 = new IGVSeekableBufferedStream(new SeekableFileStream(BAM_FILE), bufferSize);

            final int SIZE = 10000;
            final byte[] bytes1 = new byte[SIZE];
            final byte[] bytes2 = new byte[SIZE];

            reallyRead(bytes1, in1);
            reallyRead(bytes1, in1);
            in1.skip(bytes1.length);
            reallyRead(bytes1, in1);

            reallyRead(bytes2, in2);
            reallyRead(bytes2, in2);
            in2.seek(bytes2.length * 3);
            reallyRead(bytes2, in2);

            in1.close();
            in2.close();

            assertArraysEqual(bytes1, bytes2);
        }
    }

    private int reallyRead(final byte[] bytes, final IGVSeekableBufferedStream in) throws IOException {
        int read = 0, total = 0;
        do {
            read = in.read(bytes, total, bytes.length-total);
            total += read;
        } while (total != bytes.length && read > 0);

        return total;
    }


    @Test
    public void testDivisableReads()throws IOException{

        testReadsLength(1);
        testReadsLength(2);
        testReadsLength(4);
        testReadsLength(5);
        testReadsLength(10);
        testReadsLength(20);
        testReadsLength(50);
        testReadsLength(100);

    }

    private void testReadsLength(final int length) throws IOException {

        final int BUFFERED_STREAM_BUFFER_SIZE = 100;
        final byte buffer[]=new byte[BUFFERED_STREAM_BUFFER_SIZE*10];
        final SeekableFileStream fileStream = new SeekableFileStream(megabyteZerosFile);
        final IGVSeekableBufferedStream  bufferedStream = new IGVSeekableBufferedStream(fileStream,BUFFERED_STREAM_BUFFER_SIZE);

        for( int i=0; i<10*BUFFERED_STREAM_BUFFER_SIZE/length ; ++i ){
            assertEquals(bufferedStream.read(buffer, 0, length), length);
        }
    }

    private void assertArraysEqual(byte [] a1, byte [] a2) {
        assertEquals(a1.length, a2.length);
        for(int i=0; i<a1.length; i++) {
            assertEquals(a1[i], a2[i]);
        }
    }

    /**
     * Test that we can read from files larger than Integer.MAX_VALUE in length
     * @throws Exception
     */
    @Test
    public void testReadVeryLargeFile() throws Exception{
        long fileLength = ((long) Integer.MAX_VALUE) * (8);
        //long fileLength = Integer.MAX_VALUE / 2;
        long bytesToRead = fileLength;

        String path = "/dev/zero";

        SeekableStream ss = null;
        //Path doesn't exist on windows
        try{
            ss = new MaxLengthSeekableStream(path, fileLength);
        }catch(FileNotFoundException e){
            Assume.assumeNoException(e);
        }

        IGVSeekableBufferedStream bufferedStream = new IGVSeekableBufferedStream(ss);

        long bytesRead = 0;
        byte[] buf = new byte[IGVSeekableBufferedStream.DEFAULT_BUFFER_SIZE];
        while(!bufferedStream.eof() && bytesRead < bytesToRead){
            bytesRead += bufferedStream.read(buf, 0, buf.length);
            assertTrue("Could not read from source", bytesRead >= 0);
        }
        assertTrue(bytesRead >= bytesToRead);
    }



    @Test
    public void testReadStreamUnknownLength() throws Exception{
        File inFile = testFile;
        SeekableStream ss = new UnknownLengthSeekableStream(inFile);

    }

    /**
     * Class which reads from a file, but gives -1 for length.
     * Intended for checking that stream is resilient against unknown Content-Length,
     * as sometimes happens over HTTP
     */
    private static class UnknownLengthSeekableStream extends SeekableFileStream{

        UnknownLengthSeekableStream(File file) throws FileNotFoundException{
            super(file);
        }

        @Override
        public long length() {
            return -1;
        }
    }

    /**
     * Class which reads from a FileInputStream, with an artificially
     * imposed length on it.  This could act funny if {@code length} is larger
     * than the length of the actual file. Intended to be used for /dev/* classes
     * which have no actual length
     */
    private static class MaxLengthSeekableStream extends SeekableStream{

        private FileInputStream fis;
        private long position = 0;
        private long length = -1;

        /**
         *
         * @param devPath
         * @param length The /dev/* paths are unlimited, we impose a fake length
         * @throws FileNotFoundException
         */
        MaxLengthSeekableStream(String devPath, long length) throws FileNotFoundException{
            fis = new FileInputStream(devPath);
            this.length = length;
        }

        @Override
        public long length() {
            return length;
        }

        @Override
        public long position() throws IOException {
            return position;
        }

        @Override
        public void seek(long position) throws IOException {
            this.position = position;
        }

        @Override
        public int read() throws IOException {
            return fis.read();
        }

        @Override
        public int read(byte[] buffer, int offset, int length) throws IOException {
            long maxToRead = Math.min(this.length() - this.position(), (long) length);
            int iMaxToRead = Ints.saturatedCast(maxToRead);
            int bytesRead = this.fis.read(buffer, offset, iMaxToRead);
            this.position += bytesRead;
            return bytesRead;
        }

        @Override
        public void close() throws IOException {
            this.fis.close();
        }

        @Override
        public boolean eof() throws IOException {
            return position() >= length() - 1;
        }

        @Override
        public String getSource() {
            return this.fis.toString();
        }
    }

}
