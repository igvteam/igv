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

package org.broad.igv.util;

import org.junit.Test;

import java.io.File;

import static org.junit.Assert.*;

/**
 * User: jrobinso
 * Date: Mar 20, 2010
 */
public class CompressionUtilsTest {

    CompressionUtils compressionUtils = new CompressionUtils();

    @Test
    public void testCompression() {
        int sz = 1000000;
        byte[] uncompressedBytes = new byte[sz];
        for (int i = 0; i < sz; i++) {
            uncompressedBytes[i] = (byte) (Math.sin(i) * 100);
        }

        byte[] compressedBytes = compressionUtils.compress(uncompressedBytes);
        byte[] result = compressionUtils.decompress(compressedBytes);

        assertEquals(uncompressedBytes.length, result.length);
        for (int i = 0; i < result.length; i++) {
            assertEquals(uncompressedBytes[i], result[i]);
        }
    }

    @Test
    public void testCompressionChunked() {
        int sz = 1000000;
        byte[] uncompressedBytes = new byte[sz];
        for (int i = 0; i < sz; i++) {
            uncompressedBytes[i] = (byte) (Math.sin(i) * 100);
        }

        // Compress the data in 32k chunks
        int chunkSize = 32000;
        byte[] compressedBytes = compressionUtils.compress(uncompressedBytes, chunkSize);

        // Decompress.  Pass an incorrect chunk size (too small), it should still decompress correctly
        byte[] result = compressionUtils.decompress(compressedBytes, chunkSize - 1000);
        assertEquals(uncompressedBytes.length, result.length);
        for (int i = 0; i < result.length; i++) {
            assertEquals(uncompressedBytes[i], result[i]);
        }

        // Decompress.  Pass an incorrect chunk size (too large), it should still decompress correctly
        result = compressionUtils.decompress(compressedBytes, chunkSize + 1000);
        assertEquals(uncompressedBytes.length, result.length);
        for (int i = 0; i < result.length; i++) {
            assertEquals(uncompressedBytes[i], result[i]);
        }
    }

    @Test
    public void testUngzipFile_01() throws Exception {
        String inPath = TestUtils.DATA_DIR + "largegzdata.gz";
        int expSize = 23743;
        tstUngzipFile(inPath, expSize);
    }

    @Test
    public void testUngzipFile_02() throws Exception {
        String inPath = TestUtils.DATA_DIR + "testgzip.fasta.gz";
        int expSize = 775;
        tstUngzipFile(inPath, expSize);
    }

    private void tstUngzipFile(String inPath, int expSize) throws Exception {

        String outPath = inPath.substring(0, inPath.length() - 3);
        File outFile = new File(outPath);
        outFile.delete();
        outFile.deleteOnExit();
        assertFalse(outFile.exists());

        CompressionUtils.ungzipFile(inPath, outPath);
        assertTrue(outFile.exists());

        assertEquals(expSize, outFile.length());

    }

}
