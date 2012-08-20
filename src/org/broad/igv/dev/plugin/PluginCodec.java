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

package org.broad.igv.dev.plugin;

import org.broad.tribble.Feature;

import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Reads the output from a command line utility,
 * decodes it into a feature.
 * User: jacob
 * Date: 2012-Aug-01
 */
public abstract class PluginCodec implements FeatureEncoder, FeatureDecoder {

    protected List<String> cmd;
    protected Map<Argument, Object> argumentMap;
    protected Map<String, Integer> outputColumns;

    /**
     * Pattern used to split line in {@link #getNumCols(String)},
     * assuming default implementation. Subclasses may
     * simply assign a new value to columnDelimiter if counting
     * columns can be accomplished via regex
     */
    protected Pattern columnDelimiter = Pattern.compile("\\t");

    public PluginCodec(List<String> cmd, Map<Argument, Object> argumentMap) {
        this.cmd = cmd;
        this.argumentMap = argumentMap;
    }

    public void setOutputColumns(Map<String, Integer> outputColumns) {
        this.outputColumns = outputColumns;
    }

    public Feature decode(String line) {
        return decode(columnDelimiter.split(line));
    }

    public abstract Feature decode(String[] tokens);

    public int getNumCols(String line) {
        return columnDelimiter.split(line).length;
    }
}
