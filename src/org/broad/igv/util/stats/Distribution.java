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

package org.broad.igv.util.stats;

import java.io.PrintWriter;

/**
 * @author jrobinso
 * @date Jan 16, 2011
 */
public class Distribution {

    private int minBin = 0;
    private int maxBin;
    private int nBins;
    private final double[] freq;   // freq[i] = # occurences of value i
    private int n;  // # of data points

    // Create a new histogram.

    public Distribution(int maxBin) {
        this.maxBin = maxBin;
        int nBins = maxBin - minBin + 1;
        freq = new double[nBins];
    }

    public Distribution(int minBin, int maxBin) {
        this.minBin = minBin;
        this.maxBin = maxBin;
        int nBins = maxBin - minBin + 1;
        freq = new double[nBins];
    }

    // Add one occurrence of the value i.

    public void addDataPoint(int i) {
        int bin = i - minBin;
        int idx = Math.max(0, Math.min(nBins - 1, i));
        freq[idx]++;
        n++;
    }

    public double[] getDist() {
        return freq;
    }


    public void print(PrintWriter pw) {
        boolean start = false;
        for (int i = 0; i <= maxBin; i++) {
            //if (start || freq[i] > 0) {
            start = true;
            float frac = ((float) freq[i]) / n;
            pw.println((i - minBin) + "\t" + freq[i]); // + "\t" + frac);
            //}
        }
    }

}
