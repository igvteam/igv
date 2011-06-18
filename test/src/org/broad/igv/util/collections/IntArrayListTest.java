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

package org.broad.igv.util.collections;

import org.broad.igv.util.collections.IntArrayList;
import org.junit.AfterClass;

import static org.junit.Assert.assertEquals;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author jrobinso
 */
public class IntArrayListTest {
    private static final int NUM_POINTS = 10000000;

    public IntArrayListTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void test() {
        IntArrayList al = new IntArrayList(10);

        long t0 = System.currentTimeMillis();
        for (int i = 0; i < NUM_POINTS; i++) {
            al.add(i);
        }
        int[] values = al.toArray();
        assertEquals(NUM_POINTS, values.length);
        assertEquals(2, values[2]);
        System.out.println("IntArrayList time: " + (System.currentTimeMillis() - t0));
    }

    @Test
    public void testAddAll() {
        IntArrayList aList = new IntArrayList(40);
        for(int i=0; i<30; i++) {
            aList.add(i);
        }

        IntArrayList aList2 = new IntArrayList(20);
        for(int i=0; i<20; i++) {
            aList2.add(i);
        }

        aList.addAll(aList2);

        for(int i=0; i<30; i++) {
            assertEquals(i, aList.get(i));
        }
        for(int i=30; i<50; i++) {
            assertEquals(i-30, aList.get(i));
        }
    }


}