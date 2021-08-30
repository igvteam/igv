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

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.tools.FeatureSearcher;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.CancellableProgressDialog;
import org.broad.igv.ui.util.IndefiniteProgressMonitor;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.LongRunningTask;
import htsjdk.tribble.Feature;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.util.Iterator;

/**
 *
 */
class FeatureTrackUtils {

    /**
     * Attempt to load in very small chunks, we're looking for a single feature.
     */
    static int binSize = 1000;

    /**
     * Find the next/previous feature which lies outside chr:initStart-initEnd
     * @param source
     * @param chr
     * @param initStart
     * @param initEnd
     * @param forward
     * @return
     * @throws IOException
     */
    public static Feature nextFeature(FeatureSource source, String chr, int initStart, int initEnd, boolean forward) throws IOException {

        Feature f = null;

        final Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (forward) {
            // Forward
            int nextStart = initEnd;
            String nextChr = chr;
            while (nextChr != null) {
                int chrLength = genome.getChromosome(nextChr).getLength();
                while (nextStart < chrLength) {
                    int nextEnd = binSize > 0 ? nextStart + binSize : chrLength;
                    Iterator<Feature> iter = source.getFeatures(nextChr, nextStart, nextEnd);
                    if (iter != null) {
                        // The check on position should not be necessary, but not all implementations of getFeatures
                        // obey the contract to return features only in the interval.
                        while (iter.hasNext()) {
                            Feature feat = iter.next();
                            if (feat.getStart() > nextStart) {
                                return feat;
                            }
                        }
                    }
                    nextStart = nextEnd;
                }
                nextChr = genome.getNextChrName(nextChr);
                nextStart = 0;
            }
        } else {
            // Reverse
            int nextEnd = initStart;
            String nextChr = chr;
            while (nextChr != null) {
                while (nextEnd > 0) {
                    int nextStart = binSize > 0 ? Math.max(0, nextEnd - binSize) : 0;
                    Iterator<Feature> iter = source.getFeatures(nextChr, nextStart, nextEnd);
                    if (iter != null && iter.hasNext()) {
                        // The check on position should not be necessary, but not all implementations of getFeatures
                        // obey the contract to return features only in the interval.
                        Feature prevFeature = null;
                        while (iter.hasNext()) {
                            Feature feat = iter.next();
                            if (feat.getStart() < nextEnd) {
                                prevFeature = feat;
                            }
                        }
                        if (prevFeature != null) {
                            return prevFeature;
                        }
                    }
                    nextEnd = nextStart;
                }
                nextChr = genome.getPrevChrName(nextChr);
                if (nextChr != null) {
                    nextEnd = genome.getChromosome(nextChr).getLength();
                }
            }
        }

        return f;
    }

}
