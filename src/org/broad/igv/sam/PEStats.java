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

package org.broad.igv.sam;

import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;

/**
 * @author jrobinso
 * @date Mar 11, 2011
 */
public class PEStats {

    private static Logger log = Logger.getLogger(PEStats.class);

    public enum Orientation {FR, RF, F1F2, F2F1}

    ;

    private static int MAX = 10000;
    String library;
    int nPairs = 0;
    private double[] insertSizes = new double[MAX];
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
    }


    public void update(Alignment alignment) {

        if (nPairs < insertSizes.length) {
            insertSizes[nPairs] = Math.abs(alignment.getInferredInsertSize());
            nPairs++;
        }

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

        if (nPairs > 100 && insertSizes != null) {
            minThreshold = computePercentile(minPercentile);
            maxThreshold = computePercentile(maxPercentile);

            minOutlierInsertSize = computePercentile(minOutlierInsertSizePercentile);
            maxOutlierInsertSize = computePercentile(maxOutlierInsertSizePercentile);

            //log.info(library + "  " + nPairs + "  " + minThreshold + "  " + maxThreshold + " fr =" + frCount +
            //        "  ff = " + ffCount + "  rf = " + rfCount);
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
        return (int) StatUtils.percentile(insertSizes, 0, nPairs, percentile);
    }

    int getMinOutlierInsertSize() {
        return minOutlierInsertSize;
    }

    int getMaxOutlierInsertSize() {
        return maxOutlierInsertSize;
    }
}
