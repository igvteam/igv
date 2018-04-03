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

package org.broad.igv.tdf;

import junit.framework.TestCase;
import org.broad.igv.util.TestUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 * @author Jim Robinson
 * @date 5/22/12
 */
public class TDFUtilsTest extends TestCase {


    public void testTdfToBedgraph() throws Exception {

        // Create a test bedgraph file

        File testFile = new File(TestUtils.DATA_DIR + "wig/test.bedgraph");
        String trackLine = "track name=test";

        PrintWriter pw = new PrintWriter(new FileWriter(testFile));
        pw.println(trackLine);
        String chr = "chr1";
        int step = 1000;
        int span = 100;
        int start = 0;
        int end = start + step;

        int nPoints = 10000;
        float[] expectedData = new float[nPoints];
        for (int i = 0; i < nPoints; i++) {
            float expected = (float) Math.random();
            expectedData[i] = expected;
            pw.println(chr + "\t" + start + "\t" + end + "\t" + expected);
        }
        pw.close();

        // Create TDF file


        // Convert to bedgraph


        // Parse and compare results

    }
}
