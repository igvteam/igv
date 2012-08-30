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
package org.broad.igv.tools;

import org.broad.igv.track.WindowFunction;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

/**
 * @author jrobinso
 */
public class AccumulatorTest {
    protected final int numberOfPoints = (int) 1e6;
    private List<WindowFunction> wfs;
    private Map<WindowFunction, Double> values;

    public AccumulatorTest() {
    }

    @Before
    public void setUp() throws Exception {
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
        values = new HashMap<WindowFunction, Double>();
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

        // Compute stats for large number of points
        ListAccumulator accum = new ListAccumulator(wfs);
        for (int i = 0; i < numberOfPoints; i++) {
            accum.add(1, (float) Math.random());
        }
        accum.finish();

        for (WindowFunction wf : wfs) {
            double v = accum.getValue(wf);
            assertEquals(wf.getValue(), values.get(wf), v, 1.0e-2);

        }
    }


    /**
     * Pathological case,  all zeroes
     */
    @Test
    public void testZeroes() {

        // Takes too long to use all 1e6 points
        // Percentile.select is slow when all elements are the same (apparently)
        ListAccumulator accum = new ListAccumulator(wfs);
        for (int i = 0; i < 200001; i++) {
            accum.add(1, 0);
        }
        accum.finish();

        for (WindowFunction wf : wfs) {
            if (wf != WindowFunction.count) {
                double v = accum.getValue(wf);
                assertEquals(wf.getValue(), 0, v, 1.0e-2);
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
                assertTrue(wf.getValue(), Float.isNaN(v));
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
            accum.add(1, Float.NaN);
        }
        accum.finish();

        for (WindowFunction wf : wfs) {
            if (wf != WindowFunction.count) {
                float v = accum.getValue(wf);
                assertTrue(wf.getValue(), Float.isNaN(v));
            }
        }
    }

    /**
     * A single NaN
     */
    @Test
    public void testSingleNaN() {

        ListAccumulator accum = new ListAccumulator(wfs);
        accum.add(1, Float.NaN);
        accum.finish();

        for (WindowFunction wf : wfs) {
            if (wf != WindowFunction.count) {
                float v = accum.getValue(wf);
                assertTrue(wf.getValue(), Float.isNaN(v));
            }
        }
    }


    /**
     * Pathological case,  # of data points exactly equals percentile chunk size
     */
    @Test
    public void testChunkSize() {

        ListAccumulator accum = new ListAccumulator(wfs);
        for (int i = 0; i < ListAccumulator.MAX_VALUE_COUNT; i++) {
            accum.add(1, (float) Math.random());
        }
        accum.finish();
        for (WindowFunction wf : wfs) {
            double v = accum.getValue(wf);
            if (wf == WindowFunction.count) {
                assertEquals(wf.getValue(), ListAccumulator.MAX_VALUE_COUNT, v, 1.0e-2);
            } else {
                assertEquals(wf.getValue(), values.get(wf), v, 1.0e-2);
            }

        }

        accum = new ListAccumulator(wfs);
        for (int i = 0; i < ListAccumulator.MAX_VALUE_COUNT - 1; i++) {
            accum.add(1, (float) Math.random());
        }
        accum.finish();
        for (WindowFunction wf : wfs) {
            double v = accum.getValue(wf);
            if (wf == WindowFunction.count) {
                assertEquals(wf.getValue(), ListAccumulator.MAX_VALUE_COUNT - 1, v, 1.0e-2);
            } else {
                assertEquals(wf.getValue(), values.get(wf), v, 1.0e-2);
            }

        }

        accum = new ListAccumulator(wfs);
        for (int i = 0; i < ListAccumulator.MAX_VALUE_COUNT + 1; i++) {
            accum.add(1, (float) Math.random());
        }
        accum.finish();
        for (WindowFunction wf : wfs) {
            double v = accum.getValue(wf);
            if (wf == WindowFunction.count) {
                assertEquals(wf.getValue(), ListAccumulator.MAX_VALUE_COUNT + 1, v, 1.0e-2);
            } else {
                assertEquals(wf.getValue(), values.get(wf), v, 1.0e-2);
            }

        }

    }

}