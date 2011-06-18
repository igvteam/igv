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