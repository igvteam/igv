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
