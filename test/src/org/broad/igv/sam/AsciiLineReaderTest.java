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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.broad.igv.util.TestUtils;
import org.broad.tribble.readers.AsciiLineReader;
import org.junit.AfterClass;
import org.junit.BeforeClass;
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
        File testFile = new File(TestUtils.DATA_DIR + "/igv/recombRate.ens.igv.txt");


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

    @Test
    public void testSpeed() throws Exception {
        File testFile = new File(TestUtils.DATA_DIR + "/cn/HindForGISTIC.hg16.cn");


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
