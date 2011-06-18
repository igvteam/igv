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
