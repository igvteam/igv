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
 * Bin.java
 *
 * Created on October 29, 2007, 3:40 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.data;

import org.broad.igv.feature.LocusScore;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a tile bin.
 * start -- start location of the bin
 * features --  list of features that partially or fully overlap the bin
 * <p/>
 * // TODO = this is a mix of two types, a feature type bin and snp type bin.
 * //        the feature bin uses the "features" collection, snp type uses
 * startIndex / endIndex.
 * //        refactor to combine both or split.
 *
 * @author jrobinso
 */
public class FeatureBin implements Bin {

    private int start;
    private List<LocusScore> features;

    public FeatureBin(int location) {
        this.start = location;
        features = new ArrayList();
    }

    public void addFeature(LocusScore feature) {
        features.add(feature);
    }

    public List<LocusScore> getFeatures() {
        return features;
    }

    public int getFeatureCount() {

        return getFeatures().size();

    }

    /**
     * Return the scores for this bin after removing all NaN values.  If
     * there are no scores return null;
     *
     * @return
     */
    public float[] getFeatureScores() {

        int nScores = 0;
        for (LocusScore f : features) {
            if (!Float.isNaN(f.getScore())) {
                nScores++;
            }
        }

        if (nScores > 0) {
            float[] scores = new float[nScores];
            int scoreIndex = 0;
            for (LocusScore f : features) {
                float s = f.getScore();
                if (!Float.isNaN(s)) {
                    scores[scoreIndex] = s;
                    scoreIndex++;
                }
            }
            return scores;
        } else {
            return null;
        }
    }

    public int getStart() {
        return start;
    }
}

