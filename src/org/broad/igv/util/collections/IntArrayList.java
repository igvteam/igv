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

package org.broad.igv.util.collections;

/**
 * Author: jrobinso
 * Date: Jul 22, 2010
 * <p/>
 * ArrayList type collection for int types.  Purpose is to avoid the need to create an object for each entry
 * in the standard java collections.
 */
public class IntArrayList {


    private transient int[] elements;

    private int size;


    public IntArrayList() {
        this(100);
    }

    public IntArrayList(int initialCapacity) {
        if (initialCapacity < 0)
            throw new IllegalArgumentException("Illegal Capacity: " + initialCapacity);
        this.elements = new int[initialCapacity];
    }

    public IntArrayList(int[] elements) {
        this.elements = elements;
        size = elements.length;
    }

    public void add(int e) {
        if (size + 1 >= elements.length) {
            grow();
        }
        elements[size++] = e;
    }

    public void addAll(int[] args) {
        int[] newElements = new int[size + args.length];
        System.arraycopy(elements, 0, newElements, 0, size);
        System.arraycopy(args, 0, newElements, size, args.length);
        elements = newElements;
        size += args.length;
    }

    public void addAll(IntArrayList aList) {
        addAll(aList.toArray());
    }


    public int get(int idx) {
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
        int[] tmp = new int[newCapacity];
        System.arraycopy(elements, 0, tmp, 0, elements.length);
        elements = tmp;
    }


    public int[] toArray() {
        trimToSize();
        return elements;
    }


    private void trimToSize() {
        int oldCapacity = elements.length;
        if (size < oldCapacity) {
            int[] tmp = new int[size];
            System.arraycopy(elements, 0, tmp, 0, size);
            elements = tmp;
        }
    }

    public void set(int idx, int i) {
        while(idx >= elements.length) {
            grow();
        }
        elements[idx] = i;
        idx++;
        if(idx > size) size = idx;  // Tried Math.max here, it showed up in cpu profiles!
    }
}
