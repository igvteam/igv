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

package org.broad.igv.preprocess;

import org.broad.igv.data.ProcessingUtils;
import org.broad.igv.track.WindowFunction;
import org.junit.*;

import java.util.Random;

import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 */
public class ProcessingUtilsTest {

    public ProcessingUtilsTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }


    /**
     * Test of computeMedian method, of class ProcessingUtils.
     */
    @Test
    public void computePercentile() {
        Random random = new Random(94863198);
        int nPts = 10000;
        float[] values = new float[nPts];
        for (int i = 0; i < nPts; i++) {
            values[i] = random.nextFloat();
        }

        WindowFunction wf = WindowFunction.percentile10;
        float v = ProcessingUtils.computeStat(values, wf);
        assertEquals(wf.getValue(), 0.1, v, 1.0e-2);

    }
}
