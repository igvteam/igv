package org.broad.igv.charts;

/**
 * @author Jim Robinson
 * @date 10/25/11
 */
public class Axis {

    int pixelOrigin = 0;

    double min;
    double max;
    double [] ticks;
    double scale;   // Scale in value units per pixel

    int pixelWidth = 500;

    public int getPixelForValue(double value) {

           return pixelOrigin + (int) ((value - min) * scale);
    }

    public void setRange(double min, double max) {
        this.min = min;
        this.max = max;

        //NumericYTickLocator tickLocator = new NumericYTickLocator(max, min, 20, 40);
        //ticks = tickLocator.getTickMarkLocations();
        ticks = new double[] {-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
        scale =  ((double) pixelWidth) / (max - min);
    }



}
