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

import java.io.OutputStream;
import java.util.Iterator;
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
public abstract class PluginCodec<E extends Feature, D extends Feature> implements FeatureEncoder<E>, FeatureDecoder<E> {

    protected List<String> commands;

    /**
     * Pattern used to split line,
     * assuming default implementation. Subclasses may
     * simply assign a new value to columnDelimiter if counting
     * columns can be accomplished via regex
     */
    protected Pattern columnDelimiter = Pattern.compile("\\t");

    protected AsciiEncoder<E> encoder;
    protected AsciiDecoder<D> decoder;

    public PluginCodec(LineFeatureEncoder<E> lineFeatureEncoder, LineFeatureDecoder<D> lineFeatureDecoder) {
        this.encoder = new AsciiEncoder<E>(lineFeatureEncoder);
        this.decoder = new AsciiDecoder<D>(lineFeatureDecoder);
    }

    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
        this.encoder.setInputs(commands, argumentMap);
        this.decoder.setInputs(commands, argumentMap);
    }

    @Override
    public void setAttributes(List<Map<String, Object>> attributes) {
        this.decoder.setAttributes(attributes);
    }

    @Override
    public Map<String, Object> encodeAll(OutputStream outputStream, Iterator<E> features) {
        return this.encoder.encodeAll(outputStream, features);
    }

    public D decode(String line) {
        return decoder.decode(line);
    }

//    @Override
//    public String encode(E feature) {
//        return this.encoder.encode(feature);
//    }
//
//    @Override
//    public int getNumCols(String line) {
//        return this.encoder.getNumCols(line);
//    }
//
//    @Override
//    public String getHeader() {
//        return this.encoder.getHeader();
//    }
}
