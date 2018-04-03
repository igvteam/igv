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

/**
 * @author jrobinso
 */
public class StringArrayList {

    static final int maxGrowIncrement = Integer.MAX_VALUE / 10;

    int size = 0;
    String[] values;

    public StringArrayList(int maxSize) {
        values = new String[maxSize];
    }

    public void add(String v) {
        if (size >= values.length) {
            grow();
        }
        values[size] = v;
        size++;
    }

       public String get(int idx) {
        if (idx >= size) {
            throw new IndexOutOfBoundsException("Index out of bounds.  Size = " + size + " Index = " + idx);
        } else {
            return values[idx];
        }
    }

    public String[] toArray() {
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
        String[] tmp = new String[newSize];
        System.arraycopy(values, 0, tmp, 0, Math.min(tmp.length, values.length));
        values = tmp;
    }

    private void trim() {
        resize(size);
    }

}
