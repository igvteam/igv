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

import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.tribble.readers.AsciiLineReader;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Dec 6, 2009
 * Time: 1:07:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class SeekableFileStreamTest {


    @Test
    public void testSeek() throws Exception {
        String expectedLine = "chr22\t14000000\t15000000\trecombRate\t0.182272\t0.20444\t0.160104\t0\t0\t0\t0\t0\t0";
        File testFile = new File(TestUtils.DATA_DIR + "igv/recombRate.igv.txt");
        SeekableFileStream is = new SeekableFileStream(testFile);
        is.seek(149247);
        AsciiLineReader reader = new AsciiLineReader(is);
        String nextLine = reader.readLine();
        assertEquals(expectedLine, nextLine);
    }
}
