/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
