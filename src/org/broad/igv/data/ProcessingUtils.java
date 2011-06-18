/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
package org.broad.igv.data;

import org.broad.igv.util.collections.IntArrayList;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.WindowFunction;

import java.util.*;

/**
 * @author jrobinso
 */
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

    /**
     * Constructs ...
     */
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


    static class Interval {

        int start;
        int end;
        List<LocusScore> scores = new ArrayList(3);

        /**
         * Constructs ...
         *
         * @param start
         * @param end
         */
        public Interval(int start, int end) {
            this.start = start;
            this.end = end;
        }
    }

    /**
     * @param scores
     * @param wf
     * @return
     */

   
}
