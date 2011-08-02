/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.track;

import org.broad.tribble.Feature;
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
