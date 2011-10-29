package org.broad.igv.charts;

/**
 * @author Jim Robinson
 * @date 10/25/11
 */
public class Axis {

    int pixelOrigin = 0;

    double min;
    double max;
    double[] ticks;
    double scale;   // Scale in value units per pixel

    int pixelWidth = 500;

    private static final int STOP_TICK_RECURSIONS = 80;

    public int getPixelForValue(double value) {

        return pixelOrigin + (int) ((value - min) * scale);
    }

    public void setRange(double min, double max) {
        this.min = min;
        this.max = max;

        //NumericYTickLocator tickLocator = new NumericYTickLocator(max, min, 20, 40);
        //ticks = tickLocator.getTickMarkLocations();
        ticks = computeTicks(min, max, 10);
        scale = ((double) pixelWidth) / (max - min);
    }


    /**
     * Computes the tick mark start location and increment for an axis.  The start tick is the first tick to the left
     * of the minimum value.
     *
     * @param min
     * @param max
     * @param targetTickCount approximate number of ticks desired between min and max
     *
     * @return tuple -- {firstTick, tickIncrement}
     */
    public static double[] computeTicks(double min, double max, int targetTickCount) {

        if (max == min) {
            return new double[]{max};
        }

        final double v = max - min;
        double d = Math.log10(v);
        double r = Math.floor(d);

        double delta0 = Math.pow(10, r);
        int nTicks0 = (int) (v / delta0);

        double delta = delta0;
        int nTicks;

        int[] mults = {1, 2, 5, 10};
        for (int i = 0; i < mults.length; i++) {
            nTicks = nTicks0 * mults[i];
            delta = delta0 / mults[i];
            if (nTicks > 0.7 * targetTickCount) {
                break;
            }
        }

        double minTick = min - min % delta;
        if (min < 0) minTick -= delta;
        return new double[]{minTick, delta};
    }

}
