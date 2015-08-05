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

import org.apache.log4j.Logger;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.readers.PositionalBufferedStream;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * @author jacob
 * @date 2013-Jun-14
 */
public class FeatureCodecDecoder<D extends Feature> implements FeatureDecoder<D> {

    private static Logger log = Logger.getLogger(PluginSource.class);

    private FeatureCodec<D, PositionalBufferedStream> codec;

    public FeatureCodecDecoder(FeatureCodec<D, PositionalBufferedStream> codec) {
        this.codec = codec;
    }

    @Override
    public Iterator<D> decodeAll(InputStream is, boolean strictParsing) throws IOException {
        PositionalBufferedStream pis = new PositionalBufferedStream(is);

        this.codec.readHeader(pis);

        List<D> features = new ArrayList<D>();
        while (!pis.isDone()) {
            try {
                D feature = this.codec.decode(pis);
                if (feature != null) {
                    features.add(feature);
                }
            } catch (Exception e) {
                log.error(e.getMessage(), e);
                if (strictParsing) {
                    throw new RuntimeException(e);
                }
            } finally {
                pis.close();
            }
        }
        return features.iterator();
    }

    @Override
    public void setAttributes(List<Map<String, Object>> attributes) {
        //pass
    }

    @Override
    public void setInputs(List<String> commands, Map<Argument, Object> argumentMap) {
        //pass
    }

}
