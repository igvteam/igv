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


package org.broad.igv.feature.tribble;

import org.apache.log4j.Logger;
import org.broad.igv.feature.AbstractCacher;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * @author jrobinso
 * @date Jun 24, 2010
 */
public class CachingFeatureReader extends AbstractCacher implements IGVFeatureReader {

    private static Logger log = Logger.getLogger(CachingFeatureReader.class);
    private static int maxBinCount = 1000;
    private static int defaultBinSize = 16000; // <= 16 kb

    private FeatureReader tribbleFeatureReader;


    public CachingFeatureReader(FeatureReader tribbleFeatureReader) {
        this(tribbleFeatureReader, maxBinCount, defaultBinSize);
    }


    public CachingFeatureReader(FeatureReader tribbleFeatureReader, int binCount, int binSize) {
        super(binCount, binSize);
        this.tribbleFeatureReader = tribbleFeatureReader;
    }


    @Override
    protected Iterator<Feature> queryRaw(String chr, int start, int end) throws IOException {

        // Tribble iterators must be closed, so we need to copy the features and insure closure before exiting.
        CloseableTribbleIterator<Feature> iter = null;
        try {
            iter = tribbleFeatureReader.query(chr, start, end);
            List<Feature> featureList = new ArrayList<Feature>();
            while (iter.hasNext()) {
                Feature f = iter.next();
                if (f.getStart() > end) {
                    break;
                } else if (f.getEnd() < start) {
                    continue;
                } else {
                    featureList.add(f);
                }
            }
            return featureList.iterator();
        } finally {
            if (iter != null) iter.close();
        }
    }


    @Override
    public List<String> getSequenceNames() {
        return tribbleFeatureReader.getSequenceNames();
    }

    @Override
    public Object getHeader() {
        return tribbleFeatureReader.getHeader();
    }

    @Override
    public Iterator query(String chr, int start, int end) throws IOException {
        return queryCached(chr, start, end);
    }

    /**
     * Return an iterator over the entire file.  Nothing to cache,  just delegate to the wrapped reader
     *
     * @throws java.io.IOException
     */
    @Override
    public Iterator iterator() throws IOException {
        return tribbleFeatureReader.iterator();
    }


}

