package org.broad.tribble.util;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;

import static junit.framework.Assert.assertEquals;


/**
 * @author Jim Robinson
 * @date 2/18/12
 */
public class SeekableBufferedStreamTest {

    static byte[] expectedBytes;
    static File testFile;


    @BeforeClass
    public static void setup() throws Exception {
        final int fileSize = 20000;
        createTestFile(fileSize);
    }


    /**
     * Test reading individual bytes
     *
     * @throws Exception
     */
    @Test
    public void testRead() throws Exception {

        final int fileSize = expectedBytes.length;

        SeekableBufferedStream bufferedStream = new SeekableBufferedStream(new SeekableFileStream(testFile), 50);

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


    /**
     * Test reading byte arrays
     *
     * @throws Exception
     */
    @Test
    public void testReadBuffer() throws Exception {


        final int fileSize = expectedBytes.length;

        final int streamBufferSize = 50;
        SeekableBufferedStream bufferedStream = new SeekableBufferedStream(new SeekableFileStream(testFile), streamBufferSize);


        // At the end
        byte[] buffer = new byte[100];

        int pos = fileSize - buffer.length;
        bufferedStream.seek(pos);
        bufferedStream.readFully(buffer);
        for (int i = 0; i < buffer.length; i++) {
            byte b = buffer[i];
            assertEquals("" + i, expectedBytes[pos + i], b);
        }

        // Somewhere in the middle
        pos = 700;
        bufferedStream.seek(pos);
        bufferedStream.readFully(buffer);
        for (int i = 0; i < buffer.length; i++) {
            byte b = buffer[i];
            assertEquals("" + i, expectedBytes[pos + i], b);

        }

        // At the beginning
        pos = 0;
        bufferedStream.seek(pos);
        bufferedStream.readFully(buffer);
        for (int i = 0; i < buffer.length; i++) {
            byte b = buffer[i];
            assertEquals(expectedBytes[pos + i], b);
        }

        // Overlap left
        bufferedStream.seek(10000);
        bufferedStream.read();

        pos = 10000 - 75;
        bufferedStream.seek(pos);
        bufferedStream.readFully(buffer);
        for (int i = 0; i < buffer.length; i++) {
            byte b = buffer[i];
            assertEquals(expectedBytes[pos + i], b);
        }

        // Overlap right
        bufferedStream.seek(5000);
        bufferedStream.read();

        pos = 5000 + streamBufferSize - 25;
        bufferedStream.seek(pos);
        bufferedStream.readFully(buffer);
        for (int i = 0; i < buffer.length; i++) {
            byte b = buffer[i];
            assertEquals(expectedBytes[pos + i], b);
        }

        // Overlap both ends
        bufferedStream.seek(7000);
        bufferedStream.read();

        buffer = new byte[1000];
        pos = 7000 - streamBufferSize - 25;
        bufferedStream.seek(pos);
        bufferedStream.readFully(buffer);
        for (int i = 0; i < buffer.length; i++) {
            byte b = buffer[i];
            assertEquals(expectedBytes[pos + i], b);
        }

        // Completely contained
        bufferedStream.seek(3000);
        bufferedStream.read();

        buffer = new byte[20];
        pos = 3000 + 10;
        bufferedStream.seek(pos);
        bufferedStream.readFully(buffer);
        for (int i = 0; i < buffer.length; i++) {
            byte b = buffer[i];
            assertEquals(expectedBytes[pos + i], b);
        }


    }


    private static void createTestFile(int length) throws Exception {

        expectedBytes = new byte[length];

        testFile = new File("SeekableBuffered_test.bin");
        testFile.deleteOnExit();

        int range = Byte.MAX_VALUE - Byte.MIN_VALUE;
        for (int i = 0; i < length; i++) {
            expectedBytes[i] = (byte) (Byte.MIN_VALUE + (int) (Math.random() * range));
        }

        OutputStream os = null;
        try {
            os = new BufferedOutputStream(new FileOutputStream(testFile));
            os.write(expectedBytes);
        } finally {
            if (os != null) os.close();
        }

    }
}
