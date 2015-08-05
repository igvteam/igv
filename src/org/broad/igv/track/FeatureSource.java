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

package org.broad.igv.track;

import htsjdk.tribble.Feature;
import org.broad.igv.feature.LocusScore;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * A queryable source for Features.  Used by FeatureTrack
 *
 * author: jrobinso
 * Date: Jan 31, 2010
 */
public interface FeatureSource<T extends Feature> {

    /**
     * Return an iterator over all features that overlap the interval.  The coordinates are in the "UCSC" convention,
     * that is first base is zero and the interval is not inclusive of the end position,  i.e. an interval
     * spanning the first base is 0 - 1.
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    Iterator<T> getFeatures(String chr, int start, int end) throws IOException;


    /**
     * Return a list of coverage values spanning the given interval.  This can be null if coverage is not known
     * or relevant.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom   the zoom level
     * @return
     */
    List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom);


    /**
     * The featureWindowSize is the genomic interval, in bases, at which the associated feature track should start
     * rendering features.   When zoomed out beyond this interval coverage is shown, if available, or a message to
     * zoom in.
     *
     * TODO:  It seems slightly odd to have this property here, but on the other hand FeatureSource objects have
     * TODO:  access to the information required to estimate this value.  Still, consider removing this from the interface.
     *
     * @return the feature window threshold in base-pairs
     */
    int getFeatureWindowSize();

    /**
     * @param size the feature window size in base-pairs
     */

    void setFeatureWindowSize(int size);

}
