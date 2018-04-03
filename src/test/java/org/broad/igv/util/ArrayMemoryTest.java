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

package org.broad.igv.util;

import org.broad.igv.util.collections.IntArrayList;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;

import static org.junit.Assert.assertTrue;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jun 18, 2010
 * Time: 7:55:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class ArrayMemoryTest {
    private static final int size = 10000000;

    @Before
    public void setUp() throws Exception {
        System.gc();
        Runtime.getRuntime().gc();
    }

    //GetObjectSize is not particularly accurate, because it's not recursive
    //TODO 3rd party implementations of things like IntArrayList exist. Use them?
    /**NOTE:
     * RuntimeUtils has to be put in a jar and included as instrumentation in the
       junit tests, as in:
       <jvmarg value="-javaagent:${testlib.dir}/RuntimeUtils.jar"/>
      for this to work.
    **/
    @Ignore
    @Test
    public void compareMemory() throws Exception {

        IntArrayList tmp = makeIntArrayList();
        ArrayList<Integer> tmp2 = makeArrayList();
        long memIntArrList = 0;
        long memArrList = 0;

        for (int ii = 0; ii < size; ii++) {
            memArrList += JavaAgent.getObjectSize(tmp2.get(ii));
            memIntArrList += JavaAgent.getObjectSize(tmp.get(ii));
        }
        assertTrue(memIntArrList < memArrList);

    }


    public static IntArrayList makeIntArrayList() {
        IntArrayList list = new IntArrayList();
        for (int i = 0; i < size; i++) {
            list.add(i);
        }
        return list;
    }

    public static ArrayList<Integer> makeArrayList() {
        ArrayList<Integer> list = new ArrayList<Integer>();
        for (int i = 0; i < size; i++) {
            list.add(i);
        }
        return list;
    }

}
