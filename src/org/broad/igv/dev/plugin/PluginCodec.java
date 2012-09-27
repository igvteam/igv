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
 * Base class for writing/reading to/from command line.
 * Can encode/decode features.
 * <p/>
 * Type parameters are (in order):
 * encoding feature type
 * decoding feature type
 * <p/>
 * Example:
 * {@code <br>
 * public class MyClazz extends PluginCodec&lt;Feature, BasicFeature&gt;{
 * <br>
 *
 * @Override public String encode(Feature){...}
 * <br>
 * @Override public BasicFeature decode(String line){...}
 * <p/>
 * ... other methods....
 * }
 * }
 * User: jacob
 * Date: 2012-Aug-01
 */
public abstract class PluginCodec<E extends Feature, D extends Feature> implements FeatureEncoder<E>, FeatureDecoder<D> {

    protected List<String> commands;
    protected Map<Argument, Object> argumentMap;
    protected Map<String, Integer> outputColumns;

    /**
     * Pattern used to split line in {@link #getNumCols(String)},
     * assuming default implementation. Subclasses may
     * simply assign a new value to columnDelimiter if counting
     * columns can be accomplished via regex
     */
    protected Pattern columnDelimiter = Pattern.compile("\\t");

    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
        this.commands = commands;
        this.argumentMap = argumentMap;
    }

    @Override
    public void setOutputColumns(Map<String, Integer> outputColumns) {
        this.outputColumns = outputColumns;
    }

    public abstract D decode(String[] tokens);

    @Override
    public D decode(String line) {
        return decode(columnDelimiter.split(line));
    }

    @Override
    public int getNumCols(String line) {
        return columnDelimiter.split(line).length;
    }

    @Override
    public String getHeader() {
        return null;
    }

}
