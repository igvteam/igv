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

package org.broad.igv.sam;

import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.util.collections.DownsampledDoubleArrayList;

/**
 * @author jrobinso
 * @date Mar 11, 2011
 */
public class PEStats {

    private static Logger log = Logger.getLogger(PEStats.class);

    public enum Orientation {FR, RF, F1F2, F2F1}

    String library;


    //Maximum number of insertSizes to store
    private static final int MAX = 1000;
    private DownsampledDoubleArrayList insertSizes;
    private int minThreshold = 10;
    private int maxThreshold = 5000;

    // Orientation counts
    int frCount = 0;
    int rfCount = 0;

    int f1f2Count = 0;
    int f2f1Count = 0;

    Orientation orientation = null;


    /**
     * For paired arc view which are outside of the midrange
     * TODO Allow user to set?
     */
    private static final int minOutlierInsertSizePercentile = 95;
    private static final int maxOutlierInsertSizePercentile = 5;

    private int minOutlierInsertSize = minThreshold;
    private int maxOutlierInsertSize = maxThreshold;

    public PEStats(String library) {
        this.library = library;
        this.insertSizes = new DownsampledDoubleArrayList(100, MAX);
    }


    public void update(Alignment alignment) {

        insertSizes.add(Math.abs(alignment.getInferredInsertSize()));

        String po = alignment.getPairOrientation();
        if (po != null && po.length() == 4) {
            if (po.charAt(0) == 'F') {
                if (po.charAt(2) == 'F') {
                    if (po.charAt(1) == '1') {
                        f1f2Count++;
                    } else {
                        f2f1Count++;
                    }
                } else if (po.charAt(2) == 'R') {
                    frCount++;

                }
            } else if (po.charAt(0) == 'R') {
                if (po.charAt(2) == 'F') {
                    rfCount++;
                } else if (po.charAt(2) == 'R') {
                    if (po.charAt(1) == '1') {
                        f2f1Count++;
                    } else {
                        f1f2Count++;
                    }
                }
            }
        }

        // Force recomputation of orientation
        synchronized (this) {
            orientation = null;
        }
    }

    public void compute(double minPercentile, double maxPercentile) {

        if (insertSizes.size() > 100) {
            minThreshold = computePercentile(minPercentile);
            maxThreshold = computePercentile(maxPercentile);

            minOutlierInsertSize = computePercentile(minOutlierInsertSizePercentile);
            maxOutlierInsertSize = computePercentile(maxOutlierInsertSizePercentile);
        }
    }

    public int getMinThreshold() {
        return minThreshold;
    }

    public int getMaxThreshold() {
        return maxThreshold;
    }

    public synchronized Orientation getOrientation() {
        if (orientation == null) {
            int ffCount = f1f2Count + f2f1Count;
            if (ffCount > frCount && ffCount > rfCount) {
                if (f1f2Count > f2f1Count) {
                    orientation = Orientation.F1F2;
                } else {
                    orientation = Orientation.F2F1;
                }
            } else if (rfCount > frCount && rfCount > ffCount) {
                orientation = Orientation.RF;
            } else {
                orientation = Orientation.FR;
            }
        }
        return orientation;
    }

    private int computePercentile(double percentile) {
        return (int) StatUtils.percentile(insertSizes.toArray(), 0, insertSizes.size(), percentile);
    }

    int getMinOutlierInsertSize() {
        return minOutlierInsertSize;
    }

    int getMaxOutlierInsertSize() {
        return maxOutlierInsertSize;
    }
}
