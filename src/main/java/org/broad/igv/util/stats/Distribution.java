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
