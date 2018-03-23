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

import org.broad.igv.feature.BasicFeature;
import htsjdk.tribble.AsciiFeatureCodec;

/**
 * The purpose of this class is to rename incoming features, if we don't like the
 * way the external tool names them. For instance, Cufflinks uses a numbering scheme
 * This doesn't work for on-the-fly calculation since each feature will always be
 * CUFF1.1 / 1.2 etc. We would like the naming to be consistent between calculations.
 *
 * @author jacob
 * @date 2013-Apr-29
 */
public class RenameDecoder<T extends BasicFeature> extends AsciiDecoder<T>{

    public RenameDecoder(AsciiFeatureCodec<T> featureCodec){
        super(new AsciiDecoder.DecoderWrapper<T>(featureCodec));
    }

    @Override
    public T decode(String line) {
        T bf = super.decode(line);
        bf.setName(createName(bf));
        return bf;
    }

    /**
     * Create a new name for this feature.
     * Intended to be overridden, default implementation uses (1-based) start position.
     * @param bf
     * @return
     */
    protected String createName(BasicFeature bf) {
        return String.format("%d", bf.getStart() + 1);
    }
}
