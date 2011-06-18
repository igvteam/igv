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
