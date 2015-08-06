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