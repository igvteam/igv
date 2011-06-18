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
package org.broad.igv.tools;

import org.broad.igv.track.WindowFunction;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class AccumulatorTest {
    protected final int numberOfPoints = 10000000;
    List<WindowFunction> wfs;
    List<WindowFunction> percentileFunctions;
    Map<WindowFunction, Double> values;

    public AccumulatorTest() {
    }

    @Before
    public void setUp() throws Exception {
        percentileFunctions = Arrays.asList(
                WindowFunction.percentile90,
                WindowFunction.percentile10,
                WindowFunction.median,
                WindowFunction.percentile2,
                WindowFunction.percentile98);
        wfs = Arrays.asList(
                WindowFunction.min,
                WindowFunction.percentile90,
                WindowFunction.percentile10,
                WindowFunction.mean,
                WindowFunction.median,
                WindowFunction.max,
                WindowFunction.count,
                WindowFunction.percentile2,
                WindowFunction.percentile98);
        values = new HashMap();
        values.put(WindowFunction.min, 0.0);
        values.put(WindowFunction.percentile90, 0.9);
        values.put(WindowFunction.percentile10, 0.1);
        values.put(WindowFunction.mean, 0.5);
        values.put(WindowFunction.median, 0.5);
        values.put(WindowFunction.max, 1.0);
        values.put(WindowFunction.count, (double) numberOfPoints);
        values.put(WindowFunction.percentile2, 0.02);
        values.put(WindowFunction.percentile98, 0.98);


    }


    /**
     * Test calculation of all percentiles
     */
    @Test
    public void testAll() {

        // Compute stats for 1,000,000 points
        ListAccumulator accum = new ListAccumulator(wfs);
        for (int i = 0; i < numberOfPoints; i++) {
            accum.add((float) Math.random());
        }
        accum.finish();

        for (WindowFunction wf : wfs) {
            double v = accum.getValue(wf);
            assertEquals(wf.getDisplayName(), values.get(wf), v, 1.0e-2);

        }
    }


    /**
     * Pathological case,  all zeroes
     */
    @Test
    public void testZeroes() {

        // Compute stats for 1,000,000 points
        ListAccumulator accum = new ListAccumulator(wfs);
        for (int i = 0; i < numberOfPoints; i++) {
            accum.add(0);
        }
        accum.finish();

        for (WindowFunction wf : wfs) {
            if (wf != WindowFunction.count) {
                double v = accum.getValue(wf);
                assertEquals(wf.getDisplayName(), 0, v, 1.0e-2);
            }

        }
    }

    /**
     * Pathological case,  empty accumulator
     */
    @Test
    public void testEmpty() {

        ListAccumulator accum = new ListAccumulator(wfs);

        accum.finish();

        for (WindowFunction wf : wfs) {
            if (wf != WindowFunction.count) {
                float v = accum.getValue(wf);
                assertTrue(wf.getDisplayName(), Float.isNaN(v));
            }
        }
    }

    /**
     * Pathological case,  all NaN
     */
    @Test
    public void testNaN() {

        ListAccumulator accum = new ListAccumulator(wfs);
        for (int i = 0; i < 1000; i++) {
            accum.add(Float.NaN);
        }
        accum.finish();

        for (WindowFunction wf : wfs) {
            if (wf != WindowFunction.count) {
                float v = accum.getValue(wf);
                assertTrue(wf.getDisplayName(), Float.isNaN(v));
            }
        }
    }

    /**
     * A single NaN
     */
    @Test
    public void testSingleNaN() {

        ListAccumulator accum = new ListAccumulator(wfs);
        accum.add(Float.NaN);
        accum.finish();

        for (WindowFunction wf : wfs) {
            if (wf != WindowFunction.count) {
                float v = accum.getValue(wf);
                assertTrue(wf.getDisplayName(), Float.isNaN(v));
            }
        }
    }


    /**
     * Pathological case,  # of data poinst exactly equals percentile chunk size
     */
    @Test
    public void testChunkSize() {

        ListAccumulator accum = new ListAccumulator(wfs);
        for (int i = 0; i < ListAccumulator.MAX_VALUE_COUNT; i++) {
            accum.add((float) Math.random());
        }
        accum.finish();
        for (WindowFunction wf : wfs) {
            double v = accum.getValue(wf);
            if (wf == WindowFunction.count) {
                assertEquals(wf.getDisplayName(), ListAccumulator.MAX_VALUE_COUNT, v, 1.0e-2);
            } else {
                assertEquals(wf.getDisplayName(), values.get(wf), v, 1.0e-2);
            }

        }

        accum = new ListAccumulator(wfs);
        for (int i = 0; i < ListAccumulator.MAX_VALUE_COUNT - 1; i++) {
            accum.add((float) Math.random());
        }
        accum.finish();
        for (WindowFunction wf : wfs) {
            double v = accum.getValue(wf);
            if (wf == WindowFunction.count) {
                assertEquals(wf.getDisplayName(), ListAccumulator.MAX_VALUE_COUNT - 1, v, 1.0e-2);
            } else {
                assertEquals(wf.getDisplayName(), values.get(wf), v, 1.0e-2);
            }

        }

        accum = new ListAccumulator(wfs);
        for (int i = 0; i < ListAccumulator.MAX_VALUE_COUNT + 1; i++) {
            accum.add((float) Math.random());
        }
        accum.finish();
        for (WindowFunction wf : wfs) {
            double v = accum.getValue(wf);
            if (wf == WindowFunction.count) {
                assertEquals(wf.getDisplayName(), ListAccumulator.MAX_VALUE_COUNT + 1, v, 1.0e-2);
            } else {
                assertEquals(wf.getDisplayName(), values.get(wf), v, 1.0e-2);
            }

        }

    }

}