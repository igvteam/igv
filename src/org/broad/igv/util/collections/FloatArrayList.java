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

package org.broad.igv.util.collections;

/**
 * Author: jrobinso
 * Date: Jul 22, 2010
 * <p/>
 * ArrayList type collection for int types.  Purpose is to avoid the need to create an object for each entry
 * in the standard java collections.
 */
public class FloatArrayList {


    private transient float[] elements;

    private int size;


    public FloatArrayList() {
        this(100);
    }

    public FloatArrayList(int initialCapacity) {
        if (initialCapacity < 0)
            throw new IllegalArgumentException("Illegal Capacity: " + initialCapacity);
        this.elements = new float[initialCapacity];
    }

    public FloatArrayList(float[] elements) {
        this.elements = elements;
        size = elements.length;
    }

    public void add(float e) {
        if (size + 1 >= elements.length) {
            grow();
        }
        elements[size++] = e;
    }

    public void addAll(float[] args) {
        float[] newElements = new float[size + args.length];
        System.arraycopy(elements, 0, newElements, 0, size);
        System.arraycopy(args, 0, newElements, size, args.length);
        elements = newElements;
        size += args.length;
    }

    public void addAll(FloatArrayList aList) {
        addAll(aList.toArray());
    }


    public float get(int idx) {
        return elements[idx];
    }

    public int size() {
        return size;
    }

    public boolean isEmpty() {
        return size == 0;
    }

    /**
     * Empty all elements.  This logically clears the collection but does not free up any space.
     */
    public void clear() {
        size = 0;
    }

    private void grow() {
        int oldCapacity = elements.length;
        int newCapacity;
        if (oldCapacity < 10000000) {
            newCapacity = oldCapacity * 2;
        } else {
            newCapacity = (oldCapacity * 3) / 2 + 1;
        }
        float[] tmp = new float[newCapacity];
        System.arraycopy(elements, 0, tmp, 0, elements.length);
        elements = tmp;
    }


    public float[] toArray() {
        trimToSize();
        return elements;
    }


    private void trimToSize() {
        int oldCapacity = elements.length;
        if (size < oldCapacity) {
            float[] tmp = new float[size];
            System.arraycopy(elements, 0, tmp, 0, size);
            elements = tmp;
        }
    }
}
