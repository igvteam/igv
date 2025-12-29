/*
 * Represents a complete dataset (startLocations, markers, data) for a single
 * chromosome.
 */
package org.igv.data;

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
