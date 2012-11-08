/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
