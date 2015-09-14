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

import org.broad.igv.feature.Range;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Parses various cufflinks (and cuffdiff) output files as described here:
 * <p/>
 * http://cufflinks.cbcb.umd.edu/manual.html
 *
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 2:30 PM
 */
public class CufflinksParser {

    public static List<? extends Range> parse(String path) throws IOException {

        final String s = path.toLowerCase();
        if (s.endsWith("fpkm_tracking")) {
            AsciiFeatureCodec<FPKMValue> codec = new FPKMTrackingCodec(path);
            return parse(codec, path);
        } else if (s.endsWith("gene_exp.diff") || s.endsWith("cds_exp.diff")) {
            AsciiFeatureCodec<ExpDiffValue> codec = new ExpDiffCodec(path);
            return parse(codec, path);
        } else {
            throw new RuntimeException("Unsupported file type: " + path);
        }

    }

    public static <T extends Range> List<T> parse(AsciiFeatureCodec<T> codec, String path) throws IOException {
        List<T> values = new ArrayList<T>();
        AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(path, codec, false);
        Iterator<T> iter = reader.iterator();
        while(iter.hasNext()){
            values.add(iter.next());
        }
        return values;
    }

}
