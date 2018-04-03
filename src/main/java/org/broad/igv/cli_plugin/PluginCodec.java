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

package org.broad.igv.cli_plugin;

import htsjdk.tribble.Feature;

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
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap, Argument argument) {
        this.encoder.setInputs(commands, argumentMap, argument);
        this.decoder.setInputs(commands, argumentMap);
    }

    @Override
    public void setAttributes(List<Map<String, Object>> attributes) {
        this.decoder.setAttributes(attributes);
    }

    @Override
    public Map<String, Object> encodeAll(OutputStream outputStream, Iterator<? extends E> features) {
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
