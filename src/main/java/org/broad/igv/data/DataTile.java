/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

/**
 * Class is public to permit unit testing
 */
public class DataTile {

    private int[] startLocations;
    private int[] endLocations;
    private float[] values;
    private String[] featureNames;

    public DataTile(int[] startLocations, int[] endLocations, float[] values, String[] features) {
        this.startLocations = startLocations;
        this.endLocations = endLocations;
        this.values = values;
        this.featureNames = features;
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

    public float[] getValues() {
        return values;
    }

    /**
     * @return the featureNames
     */
    public String[] getFeatureNames() {
        return featureNames;
    }
}
