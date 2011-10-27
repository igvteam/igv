package org.broad.igv.charts;

import com.approximatrix.charting.coordsystem.ticklocator.NumericXTickLocator;
import com.approximatrix.charting.coordsystem.ticklocator.NumericYTickLocator;

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

        NumericYTickLocator tickLocator = new NumericYTickLocator(max, min, 20, 40);
        ticks = tickLocator.getTickMarkLocations();
        scale =  ((double) pixelWidth) / (max - min);
    }



}
