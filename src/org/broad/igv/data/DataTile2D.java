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

package org.broad.igv.data;

public class DataTile2D {

    private static DataTile2D nullDataTile;

    /**
     * Start locations
     */
    private int[] startLocations;

    /**
     * End locations -- can be null
     */
    private int[] endLocations;

    /**
     * Array of data value arrays.  Each array contains data corresponding
     * to a specific track.  So values[5] contains the data for track number 5.
     */
    private float[][] values;

    DataTile2D(int[] startLocations, int[] endLocations, float[][] values) {
        this.startLocations = startLocations;
        this.endLocations = endLocations;
        this.values = values;
    }

    public boolean isEmpty() {
        return startLocations == null || startLocations.length == 0;
    }

    public int[] getStartLocations() {
        return startLocations;
    }

    public int[] getEndLocations() {
        return endLocations;
    }

    public float[][] getValues() {
        return values;
    }

    public static DataTile2D getNullDataTile() {
        if (nullDataTile == null) {
            nullDataTile = new DataTile2D(
                    new int[]{},
                    new int[]{},
                    new float[][]{});
        }
        return nullDataTile;
    }
}
