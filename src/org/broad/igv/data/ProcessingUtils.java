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
package org.broad.igv.data;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.collections.IntArrayList;

import java.util.Arrays;
import java.util.List;

/** @deprecated Doesn't appear to be used anywhere
 * @author jrobinso
 */
@Deprecated
public class ProcessingUtils {

    protected static int[] findBoundaries(List<LocusScore> scores) {
        // Create list of all boundaries
        int[] boundaries = new int[2 * scores.size() + 1];
        for (int i = 0; i < scores.size(); i++) {
            LocusScore score = scores.get(i);
            boundaries[2 * i] = score.getStart();
            boundaries[2 * i + 1] = score.getEnd();
        }
        Arrays.sort(boundaries);
        // Remove duplicates
        IntArrayList boundaryList = new IntArrayList(boundaries.length);
        int lastPos = -1;
        for (int i = 0; i < boundaries.length; i++) {
            if (boundaries[i] != lastPos) {
                lastPos = boundaries[i];
                boundaryList.add(lastPos);
            }
        }
        boundaries = boundaryList.toArray();
        return boundaries;
    }

    private static boolean nullDataCheck(float[] data) {

        // Check for null or missing data.  Float.NaN is interpreted
        // as a placeholder for "no data.
        boolean noData = true;
        if ((data != null) && (data.length > 0)) {
            for (int i = 0; i < data.length; i++) {
                if (!Float.isNaN(data[i])) {
                    noData = false;
                    break;
                }
            }
        }


        return noData;
    }

    public ProcessingUtils() {
    }

    private static float computeQuantile(float[] data, double quantile) {

        double[] dData = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            dData[i] = data[i];
        }
        return (float) StatUtils.percentile(dData, quantile);
    }

    private static float computeMin(float[] data) {
        float min = Float.MAX_VALUE;
        for (int i = 0; i < data.length; i++) {
            if (!Float.isNaN(data[i])) {
                min = Math.min(data[i], min);
            }
        }
        return min;
    }

    private static float computeMax(float[] data) {
        float max = -Float.MAX_VALUE;
        for (int i = 0; i < data.length; i++) {
            if (!Float.isNaN(data[i])) {
                max = Math.max(data[i], max);
            }
        }
        return max;
    }

    private static float computeMean(float[] data) {
        float sum = 0.0f;
        int nPts = 0;
        for (int i = 0; i < data.length; i++) {
            if (!Float.isNaN(data[i])) {
                sum += data[i];
                nPts++;
            }
        }
        return (nPts == 0 ? Float.NaN : sum / nPts);
    }


    public static float computeStat(float[] data, WindowFunction function) {

        if (nullDataCheck(data)) {
            return Float.NaN;
        } else {
            switch (function) {
                case mean:
                    return computeMean(data);
                case median:
                    return computeQuantile(data, 50);
                case min:
                    return computeMin(data);
                case max:
                    return computeMax(data);
                case percentile2:
                    return computeQuantile(data, 2);
                case percentile10:
                    return computeQuantile(data, 10);
                case percentile90:
                    return computeQuantile(data, 90);
                case percentile98:
                    return computeQuantile(data, 98);
                case count:
                    return data.length;

            }
            return Float.NaN;
        }

    }

}
