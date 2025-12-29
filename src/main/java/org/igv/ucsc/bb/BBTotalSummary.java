package org.igv.ucsc.bb;

import org.igv.ucsc.twobit.UnsignedByteBuffer;

public class BBTotalSummary {
    public long basesCovered;
    public double minVal;
    public double maxVal;
    public double sumData;
    public double sumSquares;
    public double mean;

    static BBTotalSummary parseSummary(UnsignedByteBuffer buffer) {
        BBTotalSummary totalSummary = new BBTotalSummary();
        totalSummary.basesCovered = buffer.getLong();
        totalSummary.minVal = buffer.getDouble();
        totalSummary.maxVal = buffer.getDouble();
        totalSummary.sumData = buffer.getDouble();
        totalSummary.sumSquares = buffer.getDouble();
        totalSummary.mean = totalSummary.basesCovered > 0 ? totalSummary.sumData / totalSummary.basesCovered : 0;
        return totalSummary;
    }

    public String printString() {
        return "basesCovered " + basesCovered + "\n" +
                "minVal " + minVal + "\n" +
                "maxVal " + maxVal + "\n" +
                "sumData " + sumData + "\n" +
                "sumSquares " + sumSquares + "\n" +
                "mean " + mean + "\n";
    }
}
