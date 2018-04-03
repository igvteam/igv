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