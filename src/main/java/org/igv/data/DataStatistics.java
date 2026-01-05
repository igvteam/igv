package org.igv.data;

import org.igv.track.WindowFunction;

public class DataStatistics {

    private double mean;

    private double median;

    private double min;

    private double max;

    private double percentile10;

    private double percentile90;

    private double percentile98;

    private double stdDev;

    private boolean nullData = false;

    public static DataStatistics nullDataStat;

    static {
        nullDataStat = new DataStatistics();
        nullDataStat.nullData = true;
    }

    public double getValue(WindowFunction type) {
        if (nullData) {
            return Float.NaN;
        }
        switch (type) {
            case min:
                return min;
            case max:
                return max;
            case mean:
                return mean;
            case median:
                return median;
            case percentile10:
                return percentile10;
            case percentile90:
                return percentile90;
            case stddev:
                return stdDev;
        }
        return Float.NaN;
    }

    public double getMean() {
        return mean;
    }

    public void setMean(double mean) {
        this.mean = mean;
    }

    public double getMedian() {
        return median;
    }

    public void setMedian(double median) {
        this.median = median;
    }

    public double getMax() {
        return max;
    }

    public void setMax(double max) {
        this.max = max;
    }

    public double getPercentile90() {
        return percentile90;
    }

    public void setPercentile90(double percentile90) {
        this.percentile90 = percentile90;
    }

    public double getPercentile10() {
        return percentile10;
    }

    public void setPercentile10(double percentile10) {
        this.percentile10 = percentile10;
    }

    public double getStdDev() {
        return stdDev;
    }

    public void setStdDev(double stdDev) {
        this.stdDev = stdDev;
    }

    public double getPercentile98() {
        return percentile98;
    }

    public void setPercentile98(double percentile98) {
        this.percentile98 = percentile98;
    }

    public double getMin() {
        return min;
    }

    public void setMin(double min) {
        this.min = min;
    }
}
