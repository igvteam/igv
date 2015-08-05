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

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.ValidationStringency;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.sam.reader.WrappedIterator;

import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * @author jacob
 * @since 2012-Oct-01
 */
public class AlignmentDecoder implements FeatureDecoder<PicardAlignment> {

    @Override
    public Iterator<PicardAlignment> decodeAll(InputStream is, boolean strictParsing) throws IOException {
        SAMFileReader reader = new SAMFileReader(is);
        ValidationStringency stringency =
                strictParsing ? ValidationStringency.STRICT : ValidationStringency.SILENT;
        reader.setValidationStringency(stringency);
        return new WrappedIterator(reader.iterator());
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
