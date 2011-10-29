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
        ticks = performAutoScale(max, min, 10, 40, 0);
        scale = ((double) pixelWidth) / (max - min);
    }


    /**
     * Computes the tick mark locations on an axis
     *
     * @param max            the largest value in a data set
     * @param min            the smallest value in a data set
     * @param numTicsDesired the first guess at the number of desired tick marks
     * @return positions of tick marks along the axis
     */
    private double[] performAutoScale(double max, double min, int numTicsDesired, int maxTics, int tick_recursions) {

        // For some wiggle room on the first pass, make sure the desired number
        // of tics is at least 10
        int internalDesired = numTicsDesired;
        if (numTicsDesired < 10 && tick_recursions == 0) {
            internalDesired = 10;
        }

        double d = (max - min) / internalDesired;
        double ld = Math.log(d) / Math.log(10.0);

        // Original code
        //int ild = (int) Math.floor(ld);

        int ild = (int) Math.round(ld);

        // Axis debug output
        //System.out.println(max);
        //System.out.println(min);
        //System.out.println(numTicsDesired);

        //System.out.println(d);
        //System.out.println(ld);
        //System.out.println(ild);

        // Inrement our recursion count
        tick_recursions++;

        // Determine the number of decimal places to show
        int numDecimals = 0;
        if (ild < 0)
            numDecimals = -ild;

        double fld = Math.pow(10.0, ld - (double) ild);
        double ticValueIncrement = Math.pow(10.0, (double) ild);

        //System.out.println(ticValueIncrement);

        if (fld > 5.0) {
            ticValueIncrement *= 10.0;
            numDecimals--;
            if (numDecimals < 0)
                numDecimals = 0;
        } else if (fld > 2.0)
            ticValueIncrement = 5.0;
        else if (fld > 1.0)
            ticValueIncrement = 2.0;

        double minAdjusted = Math.floor(min / ticValueIncrement) * ticValueIncrement;
        double maxAdjusted = Math.floor(max / ticValueIncrement + 0.99999) * ticValueIncrement;
        int numTicsActual = (int) Math.floor((maxAdjusted - minAdjusted) / ticValueIncrement + 1.0e-5);

        // If simply dividing the increment by two fixes things, do it.
        if (numTicsActual > maxTics && numTicsActual / 2 <= maxTics) {
            numTicsActual = numTicsActual / 2;
            ticValueIncrement = ticValueIncrement * 2.0;
        }

        // Check for the unsolvable case here...
        if ((numTicsActual > maxTics && maxTics < 5) ||
                (tick_recursions > Math.min(STOP_TICK_RECURSIONS, numTicsDesired - 2))) {
            return simpleTics(max, min, numTicsDesired);
        } else {

            // Try again if we've exceeded the maximum number of tics
            if (numTicsActual > maxTics) {
                return performAutoScale(max, min, numTicsDesired - 1, maxTics, tick_recursions);
            }

            if (numTicsActual < numTicsDesired / 2) {
                numTicsActual = numTicsActual * 2;
                ticValueIncrement = ticValueIncrement / 2.0;
            }

            if (ticValueIncrement == 0) {
                numTicsActual = 3;
                minAdjusted = min;
                ticValueIncrement = (max - min) / 2.0;
            }

        }

        //if(numTicsActual == 1)
        //return performAutoScale(max,min,2*numTicsDesired);

        double[] return_val = new double[Math.min(numTicsActual + 1, maxTics)];

        return_val[0] = minAdjusted;
        for (int i = 1; i < return_val.length; i++) {
            return_val[i] = return_val[i - 1] + ticValueIncrement;
            //System.out.print(i);
            //System.out.print(": ");
            //System.out.println(return_val[i]);
        }
        //System.out.println();

        return return_val;

    }

    /**
     * Foolproof and simple method for determining Tic placement in cases where the tic count is
     * low theoretically.  Simply creates equidistant tic marks based on the max, min and the number
     * of tics desired inputs.  No consideration is given to how nice the tic marks will look; the
     * routine simply partitions the tic marks equally based on the span.
     *
     * @param max            the largest value to be plotted
     * @param min            the smallest value to be plotted
     * @param numTicsDesired the exact number of tics desired
     * @return unrounded positions of tick marks along the axis
     */
    private double[] simpleTics(double max, double min, int numTicsDesired) {
        if (numTicsDesired == 0) return null;
        double increment = (max - min) / ((double) (numTicsDesired - 1));
        double[] return_val = new double[numTicsDesired];

        return_val[0] = min;
        for (int i = 1; i < numTicsDesired; i++)
            return_val[i] = return_val[i - 1] + increment;

        return return_val;
    }

}
