/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.util;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * User: jrobinso
 * Date: Mar 20, 2010
 */
public class CompressionUtilsTest {


    @Test
    public void testCompression() {
        int sz = 1000000;
        byte[] uncompressedBytes = new byte[sz];
        for (int i = 0; i < sz; i++) {
            uncompressedBytes[i] = (byte) (Math.sin(i) * 100);
        }

        byte[] compressedBytes = CompressionUtils.compress(uncompressedBytes);
        byte[] result = CompressionUtils.decompress(compressedBytes);

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
        byte[] compressedBytes = CompressionUtils.compress(uncompressedBytes, chunkSize);

        // Decompress.  Pass an incorrect chunk size (too small), it should still decompress correctly
        byte[] result = CompressionUtils.decompress(compressedBytes, chunkSize - 1000);
        assertEquals(uncompressedBytes.length, result.length);
        for (int i = 0; i < result.length; i++) {
            assertEquals(uncompressedBytes[i], result[i]);
        }

        // Decompress.  Pass an incorrect chunk size (too large), it should still decompress correctly
        result = CompressionUtils.decompress(compressedBytes, chunkSize + 1000);
        assertEquals(uncompressedBytes.length, result.length);
        for (int i = 0; i < result.length; i++) {
            assertEquals(uncompressedBytes[i], result[i]);
        }
    }

}
