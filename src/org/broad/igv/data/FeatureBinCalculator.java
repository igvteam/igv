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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.data;

import org.broad.igv.feature.LocusScore;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeMap;

/**
 * @author jrobinso
 */
public class FeatureBinCalculator {

    /**
     * Allocates features to FeatureBins.  FeatureBins are equally spaced and begin at startLocation.
     *
     * @param features
     * @param nFeatureBins
     * @param startLocation
     * @param endLocation
     * @return
     */
    public List<FeatureBin> computeFeatureBins(List<? extends LocusScore> features, int nFeatureBins, double binSize, int startLocation, int endLocation) {

        TreeMap<Integer, FeatureBin> featureBins = new TreeMap();

        // Allocate features to FeatureBins
        for (LocusScore feature : features) {
            int binStart = (int) ((feature.getStart() - startLocation) / binSize);
            int binEnd = (int) ((feature.getEnd() - startLocation) / binSize);
            if (binEnd >= 0 && binStart < nFeatureBins) {
                int lastBin = Math.min(binEnd + 1, nFeatureBins);
                for (int b = Math.max(0, binStart); b < lastBin; b++) {
                    // We know this is a feature bin.  Ugly, but no time to redesign this
                    FeatureBin fBin = featureBins.get(b);
                    if (fBin == null) {
                        int location = (int) (startLocation + b * binSize);
                        fBin = new FeatureBin(location);
                        featureBins.put(b, fBin);
                    }
                    fBin.addFeature(feature);
                }
            }
        }

        Collection c = featureBins.values();
        return new ArrayList(featureBins.values());
    }


}
