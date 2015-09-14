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

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Mutation;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 *         Date: 4/9/13
 *         Time: 8:45 AM
 */
public class MutationFeatureSource implements FeatureSource<Mutation> {

    String sample;
    MutationDataManager dataManager;


    public MutationFeatureSource(String sample, MutationDataManager dataManager) {
        this.sample = sample;
        this.dataManager = dataManager;
    }

    @Override
    public Iterator<Mutation> getFeatures(String chr, int start, int end) throws IOException {
        return dataManager.getFeatures(sample, chr, start, end);
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;  //Not supported for mutation tracks
    }

    @Override
    public int getFeatureWindowSize() {
        return -1;  // Load all features for a given chromosome
    }

    @Override
    public void setFeatureWindowSize(int size) {
        // Ignored
    }


    static public class MutationDataManager {

        Range currentRange;
        Map<String, List<Mutation>> featureMap = Collections.synchronizedMap(new HashMap());
        TribbleFeatureSource tribbleFeatureSource;

        public MutationDataManager(ResourceLocator locator, Genome genome) throws IOException, TribbleIndexNotFoundException {
            this.tribbleFeatureSource = TribbleFeatureSource.getFeatureSource(locator, genome);
        }

        synchronized Iterator<Mutation> getFeatures(String trackKey, String chr, int start, int end) throws IOException {
            if (currentRange == null || !currentRange.contains(chr, start, end)) {
                Iterator<Feature> features = tribbleFeatureSource.getFeatures(chr, start, end);

                while (features.hasNext()) {
                    Mutation feat = (Mutation) features.next();
                    String thisKey = feat.getSampleId();
                    List<Mutation> keyFeatures = featureMap.get(thisKey);
                    if (keyFeatures == null) {
                        keyFeatures = new ArrayList<Mutation>();
                        featureMap.put(thisKey, keyFeatures);
                    }
                    keyFeatures.add(feat);
                    currentRange = new Range(chr, start, end);
                }

            }
            List<Mutation> featureList = featureMap.get(trackKey);
            return featureList == null ? Collections.EMPTY_LIST.iterator() : featureList.iterator();

        }
    }
}
