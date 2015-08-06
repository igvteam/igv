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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.broad.igv.util.TestUtils;
import htsjdk.tribble.readers.AsciiLineReader;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class AsciiLineReaderTest {

    public AsciiLineReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void testContent() throws Exception {
        File testFile = new File(TestUtils.DATA_DIR + "igv/recombRate.ens.igv.txt");


        AsciiLineReader reader = new AsciiLineReader(new FileInputStream(testFile));
        BufferedReader br = new BufferedReader(new FileReader(testFile));

        String arLine = null;
        int count = 0;
        while ((arLine = reader.readLine()) != null) {
            String brLine = br.readLine();
            assertEquals(arLine, brLine);
            count++;
        }
        assertTrue(count > 0);
    }

    //AsciiLineReader used to be faster than BufferedReader, no longer.
    //Still may be faster over HTTP
    //TODO Check over HTTP
    @Ignore
    @Test
    public void testSpeed() throws Exception {
        File testFile = new File(TestUtils.DATA_DIR + "cn/HindForGISTIC.hg16.cn");


        AsciiLineReader reader = new AsciiLineReader(new FileInputStream(testFile));

        long asciiCount = 0;
        long t02 = System.currentTimeMillis();
        while (reader.readLine() != null) {
            asciiCount++;
        }
        long asciiReaderTime = System.currentTimeMillis() - t02;

        BufferedReader br = new BufferedReader(new FileReader(testFile));
        long brCount = 0;
        long t0 = System.currentTimeMillis();
        while (br.readLine() != null) {
            brCount++;
        }
        long bufferedReaderTime = System.currentTimeMillis() - t0;


        // It will be considered a bug if AsciiLineReader is slower than BufferedReader
        assertTrue(bufferedReaderTime > asciiReaderTime);
        assertEquals(asciiCount, brCount);

    }
}
