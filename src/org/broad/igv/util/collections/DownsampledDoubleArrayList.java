/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
