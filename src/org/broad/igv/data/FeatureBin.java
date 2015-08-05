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

