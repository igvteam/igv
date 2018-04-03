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

package org.broad.igv.data.cufflinks;

import org.broad.igv.cli_plugin.LineFeatureDecoder;

/**
 * Wrapper class to make it easy to just extract the first sample from
 * an fpkm_tracking file.
 * @author jacob
 * @date 2013-May-28
 */
public class FPKMTrackingSampleCodec extends CufflinksCodec<FPKMSampleValue>  implements LineFeatureDecoder<FPKMSampleValue> {

    private FPKMTrackingCodec trackingCodec;

    public FPKMTrackingSampleCodec() {
        super(FPKMSampleValue.class, "Plugin");
        this.trackingCodec = new FPKMTrackingCodec(this.path);
    }

    @Override
    protected Object readHeader(String[] tokens) {
        return trackingCodec.readHeader(tokens);
    }

    @Override
    public FPKMSampleValue decode(String line) {
        FPKMValue val = trackingCodec.decode(line);
        if(val != null){
            return val.getSampleValue(0);
        }else{
            return null;
        }
    }

    @Override
    public boolean canDecode(String path) {

        String fn = path.toLowerCase();
        if(fn.endsWith(".gz")) fn = fn.substring(0, fn.length()-3);
        return fn.endsWith("fpkm_tracking");
    }
}
