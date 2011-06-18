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
 * Represents a complete dataset (startLocations, markers, data) for a single
 * chromosome.
 */
package org.broad.igv.data;

import java.util.HashMap;
import java.util.Map;

/**
 * Container for a block of data for a chromosome
 */
public class ChromosomeData {

    private String chr;

    private int[] startLocations;

    private int[] endLocations;

    private String[] probes;

    private Map<String, float[]> data;

    ChromosomeData(String chr) {
        this.chr = chr;
        data = new HashMap();
    }

    void setStartLocations(int[] locations) {
        this.startLocations = locations;
    }

    int[] getStartLocations() {
        return startLocations;
    }

    void setData(String heading, float[] x) {
        data.put(heading, x);
    }

    float[] getData(String heading) {
        return data.get(heading);
    }

    public int[] getEndLocations() {
        return endLocations;
    }

    public void setEndLocations(int[] endLocations) {
        this.endLocations = endLocations;
    }

    public String[] getProbes() {
        return probes;
    }

    public void setProbes(String[] probes) {
        this.probes = probes;
    }
}
