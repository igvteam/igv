package org.igv.util.stats;

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
