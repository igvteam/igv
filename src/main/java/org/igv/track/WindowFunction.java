/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.igv.track;

import org.igv.util.collections.CollUtils;

/**
 * @author jrobinso
 */
public enum WindowFunction implements CollUtils.Valued {

    none("None"),
    mean("Mean"),
    median("Median"),
    min("Minimum"),
    max("Maximum"),
    absoluteMax("Absolute Maximum"),
    percentile2("2nd Percentile"),
    percentile10("10th Percentile"),
    percentile90("90th Percentile"),
    percentile98("98th Percentile"),
    stddev("Standard Deviation"),
    count("Count"),
    density("Density");

    private String value = "";

    WindowFunction(String value) {
        this.value = value;
    }

    @Override
    public String getValue() {
        return value;
    }

    static public WindowFunction getWindowFunction(String name) {
        return CollUtils.findValueOf(WindowFunction.class, name);
    }

}
