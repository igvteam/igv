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

package org.broad.igv.feature;

import org.broad.igv.track.FeatureSource;
import htsjdk.tribble.Feature;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * A FeatureSource wrapper which provides caching.
 * The cache is only cleared when the source is closed.
 *
 * @author jacob
 * @date 2012/05/15
 */
public class CachingFeatureSource extends AbstractCacher implements FeatureSource {

    private static final int maxBinCount = 1000;
    private static final int defaultBinSize = 16000; // <= 16 kb

    private FeatureSource source;


    /**
     * Wraps the provided {@code source} with a caching version,
     * using default parameters.
     * @param source
     * @api
     */
    public CachingFeatureSource(FeatureSource source) {
        this(source, maxBinCount, defaultBinSize);
    }


    public CachingFeatureSource(FeatureSource source, int binCount, int binSize) {
        super(binCount, binSize);
        this.source = source;
    }

    @Override
    protected Iterator<Feature> queryRaw(String chr, int start, int end) throws IOException {
        return source.getFeatures(chr, start, end);
    }

    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        return super.queryCached(chr, start, end);
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return source.getCoverageScores(chr, start, end, zoom);
    }

    @Override
    public int getFeatureWindowSize() {
        return source.getFeatureWindowSize();
    }

    @Override
    public void setFeatureWindowSize(int size) {
        source.setFeatureWindowSize(size);
    }

    /**
     * Return the source which backs this CachingFeatureSource
     * @return
     */
    public FeatureSource getSource(){
        return this.source;
    }
}
