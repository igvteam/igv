/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.track;

import org.broad.igv.util.collections.CollUtils;

/**
 * @author jrobinso
 */
public enum WindowFunction implements CollUtils.Valued {

    none("None"),
    mean("Mean"),
    median("Median"),
    min("Minimum"),
    max("Maximum"),
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