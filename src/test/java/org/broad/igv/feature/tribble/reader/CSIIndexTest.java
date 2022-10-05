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

package org.broad.igv.feature.tribble.reader;

import org.broad.igv.util.TestUtils;
import org.junit.Before;
import org.junit.Test;

import java.io.FileInputStream;
import java.io.IOException;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 * Date: 1/23/14
 * Time: 4:24 PM
 */
public class CSIIndexTest {

    CSIIndex csi;
    CSIIndex ucsi;

    @Before
    public void init() throws IOException {


        FileInputStream fis = new FileInputStream(TestUtils.DATA_DIR + "tabix/index_test.bam.csi");
        byte[] bytes = fis.readAllBytes();
        csi = new CSIIndex();
        csi.parse(bytes, null);


        fis = new FileInputStream(TestUtils.DATA_DIR + "tabix/uncompressed_index.bam.csi");
        bytes = fis.readAllBytes();
        ucsi = new CSIIndex();
        ucsi.parse(bytes, null);


//        bai = new DiskBasedBAMFileIndex(new File(TEST_DATA_DIR, "index_test.bam.bai"), null);
//        csi = new htsjdk.samtools.CSIIndex(new File(TEST_DATA_DIR, "index_test.bam.csi"), false, null);
//        mcsi = new htsjdk.samtools.CSIIndex(new File(TEST_DATA_DIR, "index_test.bam.csi"), true, null);
//        try {
//            ucsi = new CSIIndex(Paths.get(TEST_DATA_DIR, "uncompressed_index.bam.csi"), null);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//        ubai = new DiskBasedBAMFileIndex(new File(TEST_DATA_DIR, "uncompressed_index.bam.bai"), null);
    }

    @Test
    public void testGetNumIndexLevels() {
        assertEquals(csi.getDepth(), 5);
        assertEquals(ucsi.getDepth(), 6);
    }

    @Test
    public void testGetMinShift() {
        assertEquals(csi.getMinShift(), 14);
        assertEquals(ucsi.getMinShift(), 12);
    }

    @Test
    public void testGetMaxBins() {
        assertEquals(csi.getMaxBinNumber(), 37450);
        assertEquals(ucsi.getMaxBinNumber(), 299594);
    }



}
