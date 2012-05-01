/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
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

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * User: jacob
 * Date: 2012/05/01
 */
public class CombinedFeatureSource implements org.broad.igv.track.FeatureSource {

    private FeatureSource sourceA;
    private FeatureSource sourceB;
    private Operation operation;

    private int featureWindowSize = 1000000;

    //Note: This must be the FULL path. Having bedtools on your systems path
    //is not sufficient
    static String BEDtoolsPath = "/usr/local/bin/bedtools";

    /**
     * If known, it is recommended that sourceA be the larger of the two. sourceB will
     * be loaded into memory by BEDTools.
     *
     * @param sourceA
     * @param sourceB
     * @param operation How the two sources will be combined
     */
    public CombinedFeatureSource(FeatureSource sourceA, FeatureSource sourceB, Operation operation) {
        this.sourceA = sourceA;
        this.sourceB = sourceB;
        this.operation = operation;
    }

    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        Iterator iterA = sourceA.getFeatures(chr, start, end);
        Iterator iterB = sourceB.getFeatures(chr, start, end);

        while (iterA.hasNext()) {

        }
        return null;
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null; //TODO
    }

    @Override
    public int getFeatureWindowSize() {
        return featureWindowSize;
    }

    @Override
    public void setFeatureWindowSize(int size) {
        featureWindowSize = size;
    }


    public enum Operation {
        INTERSECT("intersect"),
        MERGE("merge"),
        CLUSTER("cluster"),
        SUBTRACT("subtract");


        private String cmd;

        private Operation(String cmd) {
            this.cmd = cmd;
        }

        public String getCmd() {
            return cmd;
        }
    }
}
