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
    @Ignore
    @Test
    public void compareMemory() throws Exception {

        IntArrayList tmp = makeIntArrayList();
        ArrayList<Integer> tmp2 = makeArrayList();
        long memIntArrList = 0;
        long memArrList = 0;

        for (int ii = 0; ii < size; ii++) {
            memArrList += RuntimeUtils.getObjectSize(tmp2.get(ii));
            memIntArrList += RuntimeUtils.getObjectSize(tmp.get(ii));
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
