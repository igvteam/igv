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
package org.broad.igv.util;

public class ArrayHeapIntSorter {

    private int[] a;
    private int n;
    IntComparator c;

    public void sort(int[] a0) {
        a = a0;
        n = a.length;
        this.c = null;
        heapsort();
    }

    public void sort(int[] a0, IntComparator c) {
        a = a0;
        n = a.length;
        this.c = c;
        heapsort();
    }

    private void heapsort() {
        buildheap();
        while (n > 1) {
            n--;
            exchange(0, n);
            downheap(0);
        }
    }

    private void buildheap() {
        for (int v = n / 2 - 1; v >= 0; v--) {
            downheap(v);
        }
    }

    private void downheap(int v) {
        int w = 2 * v + 1;    // first descendant of v
        while (w < n) {
            if (w + 1 < n) // is there a second descendant?
            {
                if (c == null) {
                    if (a[w + 1] > a[w]) {
                        w++;
                    }
                } else {
                    if (c.compare(a[w + 1], a[w]) > 0) {
                        w++;
                    }
                }
            }
            // w is the descendant of v with maximum label

            if (c == null) {
                if (a[v] >= a[w]) {
                    return;
                }
            } else {
                if (c.compare(a[v], a[w]) >= 0) {
                    return;  // v has heap property
                }            // otherwise
            }
            exchange(v, w);  // exchange labels of v and w
            v = w;        // continue
            w = 2 * v + 1;
        }
    }

    private void exchange(int i, int j) {
        int t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
} 
