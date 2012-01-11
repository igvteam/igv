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