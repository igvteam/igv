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

import java.util.Random;

/**
 * Collection that holds up to maxSize doubles.  Additions > maxSize are uniformly downsampled.
 *
 * @author jrobinso
 *         Date: 7/14/12
 *         Time: 9:03 PM
 */
public class DownsampledDoubleArrayList {

    private static final Random RAND = new Random(System.currentTimeMillis());

    int maxSize;
    DoubleArrayList data;
    int downsampledCount = 0;

    /**
     * @param initialSize Initial size of the array allocated to store results
     * @param maxSize     Maximum number of doubles stored.
     */
    public DownsampledDoubleArrayList(int initialSize, int maxSize) {
        this.maxSize = maxSize;
        data = new DoubleArrayList(initialSize);
    }


    public void add(double d) {
        if (data.size() < maxSize) {
            data.add(d);
        } else {
            double samplingProb = ((double) maxSize) / (maxSize + downsampledCount + 1);
            if (RAND.nextDouble() < samplingProb) {
                int idx = (int) (RAND.nextDouble() * (data.size() - 1));
                // Replace random record with this one
                data.set(idx, d);
            }
            downsampledCount++;

        }
    }

    public double get(int idx) {
        return data.get(idx);
    }

    public int size() {
        return data.size();
    }

    public boolean isSampled() {
        return downsampledCount > 0;
    }

    public int getDownsampledCount() {
        return downsampledCount;
    }

    public double[] toArray() {
        return data.toArray();
    }
}
