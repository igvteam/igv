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

package org.broad.igv.data;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.broad.igv.util.TestUtils.DATA_DIR;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class DataUtilsTest {

    public DataUtilsTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of estimateFileMetrics method, of class DataUtils.
     */
    @Test
    public void estimateRowCount() {
        String filename = DATA_DIR + "/cn/HindForGISTIC.hg16.cn";
        int actualRowCount = 56961;

        // looking for estimate within 20$
        int estRowCount = DataUtils.estimateFileMetrics(filename).getEstRowCount();
        assertTrue(estRowCount > (int) (0.8 * actualRowCount) &&
                estRowCount < (int) (1.2 * actualRowCount));
    }

    @Test
    public void estimateProcessingTime() {
        String filename = DATA_DIR + "/cn/HindForGISTIC.hg16.cn";

        DataUtils.AsciiFileMetrics metrics = DataUtils.estimateFileMetrics(filename);

        System.out.println(DataUtils.estimatePreprocessingTime(metrics));

    }

}