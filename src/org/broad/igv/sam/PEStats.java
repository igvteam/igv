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

package org.broad.igv.sam;

import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;

import java.util.Iterator;

/**
 * @author jrobinso
 * @date Mar 11, 2011
 */
public class PEStats {

    private static Logger log = Logger.getLogger(PEStats.class);

    private static int MAX=10000;
    String rg;
    int nPairs = 0;
    private double[] insertSizes = new double[MAX];
    private int minThreshold = 10;
    private int maxThreshold = 5000;

    public PEStats(String rg) {
         this.rg = rg;
    }


    public void add(int insertSize) {
         if (nPairs < insertSizes.length) {
            insertSizes[nPairs] = Math.abs(insertSize);
            nPairs++;
        }
    }

    public void compute(double minPercentile, double maxPercentile) {

        if (nPairs > 100 && insertSizes != null) {
            minThreshold = (int) StatUtils.percentile(insertSizes, 0, nPairs, minPercentile);
            maxThreshold = (int) StatUtils.percentile(insertSizes, 0, nPairs, maxPercentile);

            log.debug(rg + "  " + nPairs + "  " + minThreshold +  "  " + maxThreshold);
        }
    }

    public int getMinThreshold() {
        return minThreshold;
    }

    public int getMaxThreshold() {
        return maxThreshold;
    }
}
