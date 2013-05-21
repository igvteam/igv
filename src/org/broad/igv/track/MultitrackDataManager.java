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

import org.broad.igv.data.Interval;
import org.broad.igv.feature.genome.Genome;
import org.broad.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * This class is designed to read a single data file to multiple tracks,
 * assuming that each track should only show data for each unique
 * value of a field. For instance, a file which had 2 samples
 * @author jacob
 */
public abstract class MultitrackDataManager<T extends Feature> {

    String path;
    Interval currentInterval;
    Map<String, List<T>> featureMap = Collections.synchronizedMap(new HashMap());
    TribbleFeatureSource tribbleFeatureSource;

    public MultitrackDataManager(String path, Genome genome) throws IOException {
        this.path = path;
        this.tribbleFeatureSource = new TribbleFeatureSource(path, genome);
    }

    synchronized Iterator<T> getFeatures(String trackKey, String chr, int start, int end) throws IOException {
        if (currentInterval == null || !currentInterval.contains(chr, start, end)) {
            Iterator<Feature> features = tribbleFeatureSource.getFeatures(chr, start, end);

            while(features.hasNext()) {
                T feat = (T) features.next();
                String thisKey = getTrackKey(feat);
                List<T> keyFeatures = featureMap.get(thisKey);
                if(keyFeatures == null) {
                    keyFeatures = new ArrayList<T>();
                    featureMap.put(thisKey, keyFeatures);
                }
                keyFeatures.add(feat);
                currentInterval = new Interval(chr, start, end);
            }

        }
        List<T> featureList = featureMap.get(trackKey);
        return featureList == null ? Collections.EMPTY_LIST.iterator() : featureList.iterator();

    }

    /**
     * Retrieve the from this feature which determines which track it belongs to.
     * e.g. mutation.getSampleId()
     * @param feat
     * @return
     */
    protected abstract String getTrackKey(T feat);

}
