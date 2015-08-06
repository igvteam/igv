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

package org.broad.igv.sam.reader;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.Alignment;
import org.broad.igv.tools.IGVToolsTest;
import org.broad.igv.util.TestUtils;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import java.io.File;
import java.util.Iterator;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012/03/01
 */
public class AlignmentReaderFactoryTest extends AbstractHeadlessTest {

    @Rule
    public TestRule testTimeout = new Timeout((int) 1e3 * 60);

    @Test
    public void testGetMergedReader() throws Exception {
        String relfiname = "HG00171.hg18.bam";
        String path = TestUtils.LARGE_DATA_DIR + "temprel.bam.list";
        String[] relnames = IGVToolsTest.generateRepLargebamsList(path, relfiname, 3, false);
        for (String s : relnames) {
            assertFalse((new File(s)).isAbsolute());
        }

        String apath = TestUtils.LARGE_DATA_DIR + "tempabs.bam.list";
        String[] absnames = IGVToolsTest.generateRepLargebamsList(apath, relfiname, 3, true);
        for (String s : absnames) {
            assertTrue((new File(s)).isAbsolute());
        }

        AlignmentReader relreader = AlignmentReaderFactory.getReader(path, false);
        AlignmentReader absreader = AlignmentReaderFactory.getReader(apath, false);
        Iterator<Alignment> iter = relreader.iterator();
        Iterator<Alignment> absiter = absreader.iterator();

        while (iter.hasNext()) {
            Alignment relA = iter.next();
            Alignment absA = absiter.next();
            //Only check 1 field for speed, these are long strings so unlikely to be equal by accident
            assertEquals(absA.getCigarString(), relA.getCigarString());
        }
    }

    public static void main(String[] args) {
        String path = "http://www.example.com";
        File fi = new File(path);
        System.out.println(fi.isAbsolute());
    }


}
