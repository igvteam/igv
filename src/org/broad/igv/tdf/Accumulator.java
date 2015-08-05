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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.DownsampledDoubleArrayList;

import java.util.HashSet;
import java.util.Set;

/**
 * Accumulator for a single window function.
 * <p/>
 * Estimating percentiles -- weighted average of multiple estimates
 *
 * @author jrobinso
 */
public class Accumulator {

    private static Logger log = Logger.getLogger(Accumulator.class);

    private static int MAX_VALUE_COUNT = 100000;

    boolean isFinished = false;
    WindowFunction windowFunction;
    float sum = 0.0f;
    int basesCovered = 0;
    int nPts = 0;
    float value = Float.NaN;

    DownsampledDoubleArrayList valueList;  // List used to accumulate values for percentile calculations


    // Optional -- keep some representative data and probe names for popup text
    int nRepValues;
    float[] repData;
    String[] repProbes;

    public Accumulator(WindowFunction windowFunction, int nRepValues) {
        this(windowFunction);
        if (nRepValues > 0) {
            this.nRepValues = nRepValues;
            this.repData = new float[nRepValues];
            this.repProbes = new String[nRepValues];
        }
    }


    public Accumulator(WindowFunction windowFunction) {
        this.windowFunction = windowFunction;
        if (PERCENTILE_WINDOW_FUNCTIONS.contains(windowFunction)) {
            valueList = new DownsampledDoubleArrayList(100, MAX_VALUE_COUNT);
        }
    }

    public boolean hasData() {
        return basesCovered > 0;
    }

    public void add(int nBases, float v, String probe) {

        if (isFinished) {
            log.error("Attempt to add data to a finalized accumulator");
            throw new RuntimeException("Attempt to add data to a finalized accumulator");
        }

        // Some older TDF files created in previous versions of igvtools from wig files were improperly coded,
        // with start=end, resulting in an nBases value of zero.  This is not a possible value, so threshold it
        if (nBases < 1) nBases = 1;

        if (!Float.isNaN(v)) {
            if (repData != null && nPts < repData.length) {
                repData[nPts] = v;
                repProbes[nPts] = probe;
            }
            switch (windowFunction) {
                case min:
                    value = Float.isNaN(value) ? v : Math.min(value, v);
                    break;
                case max:
                    value = Float.isNaN(value) ? v : Math.max(value, v);
                    break;
                case mean:
                    sum += nBases * v;
                    break;
                default:
                    if (valueList != null) {
                        valueList.add(v);
                    }
            }
            nPts++;
            basesCovered += nBases;
        }
    }


    public void finish() {

        if (isFinished) {
            return;
        }

        if (windowFunction == WindowFunction.mean) {
            value = Float.isNaN(sum) ? Float.NaN : sum / basesCovered;
        } else if (valueList != null) {
            if (valueList.size() == 0) {
                value = Float.NaN;
            } else if (valueList.size() == 1) {
                value = (float) valueList.get(0);
            } else {
                double[] valueArray = valueList.toArray();
                double p = this.getPercentile(windowFunction);
                if (p > 0) {
                    value = (float) StatUtils.percentile(valueArray, p);
                } else {
                    value = Float.NaN;
                }
            }

        }

        valueList = null;
        isFinished = true;

    }


    public int getNpts() {
        return nPts;
    }

    public String[] getRepProbes() {
        return repProbes;
    }

    public float getValue() {
        if (!isFinished) finish();
        return value;
    }


    public double getPercentile(WindowFunction wf) {
        switch (wf) {
            case percentile2:
                return 2;
            case percentile10:
                return 10;
            case percentile90:
                return 90;
            case percentile98:
                return 98;
            case median:
                return 50;
            default:
                return -1.0;
        }
    }

    public float[] getRepData() {
        return repData;
    }


    class PercentileValue {
        int nPoints;
        double value;

        PercentileValue(int nPoints, double value) {
            this.nPoints = nPoints;
            this.value = value;
        }
    }

    static Set<WindowFunction> PERCENTILE_WINDOW_FUNCTIONS = new HashSet();

    static {
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.median);
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.percentile2);
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.percentile10);
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.percentile90);
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.percentile98);
    }

}
