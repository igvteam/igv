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

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

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
}
