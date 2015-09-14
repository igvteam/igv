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

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.util.ParsingUtils;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.LineIterator;

/**
 * @author jacob
 * @date 2013-Apr-18
 */
public abstract class CufflinksCodec<T extends Feature> extends AsciiFeatureCodec<T> {

    private static Logger log = Logger.getLogger(CufflinksCodec.class);

    String path;

    protected CufflinksCodec(Class<T> clazz, String path){
        super(clazz);
        this.path = path;
    }

    protected abstract Object readHeader(String[] tokens);

    @Override
    public Object readActualHeader(LineIterator reader){
        String headerLine = null;
        try {
            headerLine = reader.next();
            String[] tokens = ParsingUtils.TAB_PATTERN.split(headerLine);
            return readHeader(tokens);
        } catch (Exception e) {
            log.error(e.getMessage(), e);
            throw new DataLoadException("Error reading header: " + e.getMessage(), this.path);
        }
    }
}

