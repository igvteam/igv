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

package org.broad.igv.sam.reader;


import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.peaks.Peak;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.FeatureWrappedAlignment;
import org.broad.igv.feature.tribble.PSLCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;

/**
 * @author jrobinso
 * @date Aug 5, 2010
 */

public class PSLAlignmentParser implements AlignmentParser {

    PSLCodec codec;

    public PSLAlignmentParser() {
        codec = new PSLCodec(); //genome);
    }

    public Alignment readNextRecord(AsciiLineReader reader) throws IOException {

        String nextLine;
        while ((nextLine = reader.readLine()) != null) {

            Feature f = codec.decode(nextLine);
            if(f != null && f instanceof BasicFeature) {
                return new FeatureWrappedAlignment((BasicFeature) f);
            }
        }
        return null;

    }
}
