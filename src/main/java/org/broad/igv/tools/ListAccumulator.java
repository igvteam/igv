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
package org.broad.igv.tools;

import org.broad.igv.util.collections.DoubleArrayList;
import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;
import org.broad.igv.track.WindowFunction;

import java.util.*;

/**
 * Estimating percentiles -- weighted average of multiple estimates
 *
 * @author jrobinso
 */
public class ListAccumulator {

    static Set<WindowFunction> PERCENTILE_WINDOW_FUNCTIONS = new HashSet();
    public static int MAX_VALUE_COUNT = 100000;
    private static Logger log = Logger.getLogger(ListAccumulator.class);

    static {
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.median);
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.percentile2);
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.percentile10);
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.percentile90);
        PERCENTILE_WINDOW_FUNCTIONS.add(WindowFunction.percentile98);
    }


    boolean isFinished = false;

    List<WindowFunction> windowFunctions;
    List<WindowFunction> quantileFunctions;
    Map<WindowFunction, List<PercentileValue>> percentiles = new HashMap();
    DoubleArrayList values = null;
    float sum = 0.0f;
    int basesCovered = 0;
    int nPts = 0;

    float min = Float.NaN;
    float max = Float.NaN;
    float mean = Float.NaN;
    float median = Float.NaN;
    float percentile2 = Float.NaN;
    float percentile10 = Float.NaN;
    float percentile90 = Float.NaN;
    float percentile98 = Float.NaN;


    public ListAccumulator(Collection<WindowFunction> windowFunctions) {
        this.windowFunctions = new ArrayList(windowFunctions);
        quantileFunctions = new ArrayList();
        for (WindowFunction wf : windowFunctions) {
            if (PERCENTILE_WINDOW_FUNCTIONS.contains(wf)) {
                quantileFunctions.add(wf);
                if (values == null) {
                    values = new DoubleArrayList();
                }
            }
        }
    }

    public void add(int w, float v) {
        if (!Float.isNaN(v)) {
            min = Float.isNaN(min) ? v : Math.min(min, v);
            max = Float.isNaN(max) ? v : Math.max(max, v);
            sum += w*v;
            basesCovered +=w;
            nPts++;
            if (values != null) {
                values.add(v);
                if (values.size() > MAX_VALUE_COUNT) {
                    computePercentiles();
                    values.clear();
                }
            }
        }
    }


    public void finish() {

        if (isFinished) {
            return;
        }

        mean = Float.isNaN(sum) ? Float.NaN : sum / basesCovered;

        if (values != null) {
            if (nPts == 1) {
                for (WindowFunction wf : quantileFunctions) {
                    setValue(wf, mean);
                }
            } else {
                if (values.size() > 1) {
                    computePercentiles();
                }
                for (WindowFunction wf : quantileFunctions) {

                    List<PercentileValue> pList = percentiles.get(wf);
                    float v = Float.NaN; // <= Default,
                    if (pList != null && pList.size() > 0) {
                        double weightedSum = 0;
                        double sumOfWeights = 0;
                        for (PercentileValue pv : pList) {
                            double weight = (double) pv.nPoints / nPts;
                            sumOfWeights += weight;
                            weightedSum += weight * pv.value;
                        }
                        v = (float) (weightedSum / sumOfWeights);
                    }
                    setValue(wf, v);

                }

            }
        }
        values = null;
        isFinished = true;

    }

    private void computePercentiles() {
        if (values != null) {
            double[] valueArray = values.toArray();
            for (WindowFunction wf : quantileFunctions) {
                double p = this.getPercentile(wf);
                if (p > 0) {
                    float v = (float) StatUtils.percentile(valueArray, p);
                    if (Float.isInfinite(v)) {
                        log.error("Infinite percentile (" + wf + ")");
                    } else {
                        List<PercentileValue> pList = percentiles.get(wf);
                        if (pList == null) {
                            pList = new ArrayList();
                            percentiles.put(wf, pList);
                        }
                        pList.add(new PercentileValue(valueArray.length, v));
                    }
                }
            }
        }

    }

    private void setValue(WindowFunction wf, float value) {
        switch (wf) {
            case mean:
                mean = value;
                break;
            case median:
                median = value;
                break;
            case min:
                min = value;
                break;
            case max:
                max = value;
                break;
            case percentile2:
                percentile2 = value;
                break;
            case percentile10:
                percentile10 = value;
                break;
            case percentile90:
                percentile90 = value;
                break;
            case percentile98:
                percentile98 = value;
                break;
            default:
                System.err.println("Unexpected window function: " + wf.toString());
        }


    }

    public float getValue(WindowFunction wf) {

        switch (wf) {
            case mean:
                return mean;
            case median:
                return median;
            case min:
                return min;
            case max:
                return max;
            case percentile2:
                return percentile2;
            case percentile10:
                return percentile10;
            case percentile90:
                return percentile90;
            case percentile98:
                return percentile98;
            case count:
                return nPts;
            default:
                System.err.println("Unexpected window function: " + wf.toString());

        }
        return Float.NaN;

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

}
