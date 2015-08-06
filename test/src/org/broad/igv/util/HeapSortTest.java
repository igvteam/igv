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
package org.broad.igv.util;

import org.junit.AfterClass;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.Comparator;

/**
 * @author jrobinso
 */
public class HeapSortTest {

    public HeapSortTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void tesArraySort() {

        int nPts = 1000000;

        int[] array = new int[nPts];

        IntComparator comp = new IntComparator() {

            public int compare(int arg0, int arg1) {
                return arg0 - arg1;
            }
        };

        //////////////////////////////////////
        for (int i = 0; i < nPts; i++) {
            array[i] = (int) (nPts * Math.random());
        }
        long t0 = System.currentTimeMillis();
        ArrayHeapIntSorter sorter = new ArrayHeapIntSorter();
        sorter.sort(array, comp);
        System.out.println("Random sort: " + (System.currentTimeMillis() - t0));
        assertOrder(array);


        for (int i = 0; i < nPts; i++) {
            array[i] = (int) (nPts * Math.random());
        }
        t0 = System.currentTimeMillis();
        sorter = new ArrayHeapIntSorter();
        sorter.sort(array, new IntComparator() {

            public int compare(int arg0, int arg1) {
                return arg1 - arg0;
            }
        });
        System.out.println("Random sort (reversed): " + (System.currentTimeMillis() - t0));
        assertOrderReversed(array);


        Integer[] objArray = new Integer[nPts];
        for (int i = 0; i < nPts; i++) {
            objArray[i] = (int) (nPts * Math.random());
        }
        t0 = System.currentTimeMillis();
        Arrays.sort(objArray, new Comparator<Integer>() {

            public int compare(Integer arg0, Integer arg1) {
                return arg0.intValue() - arg1.intValue();
            }
        });
        System.out.println("Random sort (Arrays): " + (System.currentTimeMillis() - t0));


        for (int i = 0; i < nPts; i++) {
            array[i] = i;
        }
        t0 = System.currentTimeMillis();

        sorter.sort(array, comp);
        System.out.println("Sorted sort: " + (System.currentTimeMillis() - t0));


        for (int i = 0; i < nPts; i++) {
            array[i] = nPts - i;
        }

        t0 = System.currentTimeMillis();
        sorter.sort(array, comp);
        System.out.println("Reverse sorted sort: " + (System.currentTimeMillis() - t0));

        // Almost sorted
        for (int i = 0; i < nPts; i++) {
            if (i % 10 == 0) {
                array[i] = (int) (nPts * Math.random());
            } else {
                array[i] = i;
            }
        }
        t0 = System.currentTimeMillis();
        sorter.sort(array, comp);
        System.out.println("Almost sorted sort: " + (System.currentTimeMillis() - t0));
        assertOrder(array);
    }

    private static void assertOrder(int[] array) {
        int last = Integer.MIN_VALUE;
        for (int i : array) {
            assertTrue("last=" + last + " i=" + i, i >= last);
            last = i;
        }
    }

    private static void assertOrderReversed(int[] array) {
        int last = Integer.MAX_VALUE;
        for (int i : array) {
            assertTrue("last=" + last + " i=" + i, i <= last);
            last = i;
        }
    }
}