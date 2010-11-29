/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jun 18, 2010
 * Time: 7:55:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class ArrayMemoryTest {


    public static void main(String[] args) throws IOException {

        char c;

        do {
            IntArrayList tmp = makeList();
            Runtime.getRuntime().gc();
            System.out.println(RuntimeUtils.getAvailableMemory());
            c = (char) System.in.read();

        }
        while (c != 'x');


    }


    public static IntArrayList makeList() {
        IntArrayList list = new IntArrayList();
        for (int i = 0; i < 10000000; i++) {
            list.add(i);
        }
        return list;
    }
}
