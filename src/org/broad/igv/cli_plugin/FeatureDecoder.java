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

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * User: jacob
 * Date: 2012-Aug-02
 *
 * @see PluginCodec
 * @see FeatureEncoder
 */
public interface FeatureDecoder<T extends Feature> {


    /**
     * Decode all features from the specified input stream
     *
     * @param is
     * @param strictParsing If true, errors are thrown if we cannot parse a given line.
     *                      If false, we simply skip that line. Mostly applicable to AsciiDecoders,
     *                      but could be useful more generally.
     * @return Iterator of decoded features
     * @throws IOException
     */
    Iterator<T> decodeAll(InputStream is, boolean strictParsing) throws IOException;

    /**
     * @param attributes List of maps containing attributes which were provided by
     *                   {@link FeatureEncoder#encodeAll(java.io.OutputStream, java.util.Iterator)}
     */
    void setAttributes(List<Map<String, Object>> attributes);


    /**
     * It may be the case that the output is processed differently, depending on the input.
     * We allow for that by providing the commands and arguments here. Implementations are not
     * required to do anything with the input
     *
     * @param commands    Command portions of the input, e.g. {"find", "."}
     * @param argumentMap Arguments with their values. e.g. "-name", "myFile"
     */
    void setInputs(List<String> commands, Map<Argument, Object> argumentMap);

}
