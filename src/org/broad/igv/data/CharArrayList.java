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

package org.broad.igv.data;

/**
 * @author jrobinso
 */
public class CharArrayList {

    static final int maxGrowIncrement = Integer.MAX_VALUE / 10;

    int size = 0;
    char[] values;

    public CharArrayList(int maxSize) {
        values = new char[maxSize];
    }

    public void add(char v) {
        if (size >= values.length) {
            grow();
        }
        values[size] = v;
        size++;
    }

    public char[] toArray() {
        trim();
        return values;
    }

    private void grow() {
        if (values.length >= Integer.MAX_VALUE) {
            throw new RuntimeException("Maximum array size exceeded");
        }
        int increment = (int) (Math.max(1000, 0.2 * values.length));
        int newSize = Math.min(Integer.MAX_VALUE, values.length + increment);
        resize(newSize);

    }

    private void resize(int newSize) {
        char[] tmp = new char[newSize];
        System.arraycopy(values, 0, tmp, 0, Math.min(tmp.length, values.length));
        values = tmp;
    }

    private void trim() {
        resize(size);
    }

}
