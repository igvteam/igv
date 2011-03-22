/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import org.broad.igv.util.collections.DoubleArrayList;
import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.track.WindowFunction;

import java.util.*;

/**
 * Accumulator for a single window function.
 * <p/>
 * Estimating percentiles -- weighted average of multiple estimates
 *
 * @author jrobinso
 */
public class Accumulator {

    public static int MAX_VALUE_COUNT = 100000;
    private static Logger log = Logger.getLogger(Accumulator.class);

    boolean isFinished = false;

    WindowFunction windowFunction;
    List<PercentileValue> percentiles;
    float sum = 0.0f;
    int nPts = 0;
    float value = Float.NaN;
    DoubleArrayList values;


    public Accumulator(WindowFunction windowFunction) {
        this.windowFunction = windowFunction;
        if (PERCENTILE_WINDOW_FUNCTIONS.contains(windowFunction)) {
            values = new DoubleArrayList();
            percentiles = new ArrayList();
        }
    }

    public void add(float v) {
        if (!Float.isNaN(v)) {
            switch (windowFunction) {
                case min:
                    value = Float.isNaN(value) ? v : Math.min(value, v);
                    break;
                case max:
                    value = Float.isNaN(value) ? v : Math.max(value, v);
                    break;
                case mean:
                    sum += v;
                    break;
                default:
                    if (values != null) {
                        values.add(v);
                        if (values.size() > MAX_VALUE_COUNT) {
                            computePercentiles();
                            values.clear();
                        }
                    }
            }
        }
        nPts++;
    }


    public void finish() {

        if (isFinished) {
            return;
        }

        if (windowFunction == WindowFunction.mean) {
            value = Float.isNaN(sum) ? Float.NaN : sum / nPts;
        } else if (values != null) {
            if (values.size() == 1) {
                value = (float) values.get(0);
            } else {
                if (values.size() > 1) {
                    computePercentiles();
                }


                float v = Float.NaN; // <= Default,
                if (percentiles != null && percentiles.size() > 0) {
                    double weightedSum = 0;
                    double sumOfWeights = 0;
                    for (PercentileValue pv : percentiles) {
                        double weight = (double) pv.nPoints / nPts;
                        sumOfWeights += weight;
                        weightedSum += weight * pv.value;
                    }
                    v = (float) (weightedSum / sumOfWeights);
                }
                value = v;


            }
        }
        values = null;
        isFinished = true;

    }

    private void computePercentiles() {
        if (values != null) {
            double[] valueArray = values.toArray();
            double p = this.getPercentile(windowFunction);
            if (p > 0) {
                float v = (float) StatUtils.percentile(valueArray, p);
                if (Float.isInfinite(v)) {
                    log.error("Infinite percentile (" + windowFunction + ")");
                } else {
                    percentiles.add(new PercentileValue(valueArray.length, v));
                }
            }
        }


    }


    public float getValue() {
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
