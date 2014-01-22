/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.Mutation;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;

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
