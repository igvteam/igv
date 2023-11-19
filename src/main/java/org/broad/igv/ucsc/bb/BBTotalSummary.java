package org.broad.igv.ucsc.bb;

import org.broad.igv.ucsc.UnsignedByteBuffer;

public class BBTotalSummary {
    public long basesCovered;
    public double minVal;
    public double maxVal;
    public double sumData;
    public double sumSquares;
    public double mean;
    public double stddev;

    static BBTotalSummary parseSummary(UnsignedByteBuffer buffer) {
        BBTotalSummary totalSummary = new BBTotalSummary();
        totalSummary.basesCovered = buffer.getLong();
        totalSummary.minVal = buffer.getDouble();
        totalSummary.maxVal = buffer.getDouble();
        totalSummary.sumData = buffer.getDouble();
        totalSummary.sumSquares = buffer.getDouble();
        totalSummary.computeStats();
        return totalSummary;
    }

    void computeStats() {

        long n = this.basesCovered;
        if (n > 0) {
            this.mean = this.sumData / n;
            this.stddev = Math.sqrt(this.sumSquares / (n - 1));
        }
    }
}
