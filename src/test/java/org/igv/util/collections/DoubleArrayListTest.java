/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.util.collections;

import org.igv.util.collections.IntArrayList;
import org.junit.AfterClass;

import static org.junit.Assert.assertEquals;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author jrobinso
 */
public class DoubleArrayListTest {
    private static final int NUM_POINTS = 10000;

    public DoubleArrayListTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void test() {
        DoubleArrayList al = new DoubleArrayList(10);

        long t0 = System.currentTimeMillis();
        for (int i = 0; i < NUM_POINTS; i++) {
            al.add(i);
        }
        double[] values = al.toArray();
        assertEquals(NUM_POINTS, values.length);
        assertEquals(2, values[2], 1.0e-6);
        System.out.println("IntArrayList time: " + (System.currentTimeMillis() - t0));
    }

    @Test
    public void testAddAll() {
        DoubleArrayList aList = new DoubleArrayList(40);
        for(int i=0; i<30; i++) {
            aList.add(i);
        }

        DoubleArrayList aList2 = new DoubleArrayList(20);
        for(int i=0; i<20; i++) {
            aList2.add(i);
        }

        aList.addAll(aList2);

        for(int i=0; i<30; i++) {
            assertEquals(i, aList.get(i), 1.0e-6);
        }
        for(int i=30; i<50; i++) {
            assertEquals(i-30, aList.get(i), 1.0e-6);
        }
    }


}