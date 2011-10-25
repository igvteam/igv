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
 * ExpressionDataset.java
 *
 * Created on October 18, 2007, 2:20 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.data.expression;

import org.broad.igv.Globals;
import org.broad.igv.data.Dataset;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;

import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class ExpressionDataset implements Dataset {

    private String name;
    private TrackType type = TrackType.GENE_EXPRESSION;
    private Genome genome;
    private String[] columnHeadings;
    private boolean normalized = false;
    private boolean logValues = false;

    /**
     * Map colum heading -> index for effecient reverse lookup
     */
    private Map<String, Integer> headingIndexMap = new HashMap();

    Map<String, int[]> startLocationMap = new HashMap();
    Map<String, int[]> endLocationMap = new HashMap();

    private Map<String, Integer> longestFeatureMap;

    /**
     * Map of chromosome -> array of data values
     */
    Map<String, Map<String, float[]>> dataMap = new HashMap();

    /**
     * Map of chromosome -> array of feature names
     */
    Map<String, String[]> featureNameMap = new HashMap();

    private TrackProperties trackProperties = new TrackProperties();

    /**
     * Creates a new instance of ExpressionDataset
     */
    public ExpressionDataset(Genome genome) {
        this.genome = genome;
    }


    // Todo -- implement

    public float getDataMin() {
        return -3f;
    }

    public float getDataMax() {
        return 3f;
    }


    public void setColumnHeadings(String[] columnHeadings) {
        this.columnHeadings = columnHeadings;
        for (int i = 0; i < columnHeadings.length; i++) {
            headingIndexMap.put(columnHeadings[i], i);
        }
    }

    public String[] getTrackNames() {
        return columnHeadings;
    }

    public String getName() {
        return name;
    }

    public void setType(TrackType type) {
        this.type = type;
    }


    public TrackType getType() {
        return type;
    }

    public boolean isEmpty() {
        return startLocationMap.isEmpty();
    }


    public String[] getChromosomes() {
        return startLocationMap.keySet().toArray(new String[0]);
    }


    public void setFeatureNames(String chr, String[] names) {
        this.featureNameMap.put(chr, names);
    }

    public String[] getFeatureNames(String chr) {
        return featureNameMap.get(chr);
    }

    public void setStartLocations(String chr, int[] startLocations) {
        this.startLocationMap.put(chr, startLocations);
    }

    public int[] getStartLocations(String chr) {
        return startLocationMap.get(chr);
    }

    public void setEndLocations(String chr, int[] endLocations) {
        this.endLocationMap.put(chr, endLocations);
    }

    public int[] getEndLocations(String chr) {

        if (chr.equals(Globals.CHR_ALL)) {
            return null;
        }

        return endLocationMap.get(chr);
    }

    public boolean isLogNormalized() {
        return normalized;
    }

    public void setData(String heading, String chr, float[] data) {
        Map<String, float[]> tmp = dataMap.get(heading);
        if (tmp == null) {
            tmp = new HashMap();
            dataMap.put(heading, tmp);
        }
        tmp.put(chr, data);
    }

    public float[] getData(String heading, String chr) {
        Map<String, float[]> tmp = dataMap.get(heading);
        if (tmp != null) {
            return tmp.get(chr);
        }
        return null;
    }

    public void setName(String name) {
        this.name = name;
    }

    public boolean isLogValues() {
        return logValues;
    }

    public void setLogValues(boolean logValues) {
        this.logValues = logValues;
    }

    public void setNormalized(boolean normalized) {
        this.normalized = normalized;
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    public void setTrackLine(String trackLine) {
        ParsingUtils.parseTrackLine(trackLine, trackProperties);
    }

    public Integer getLongestFeature(String chr) {
        return longestFeatureMap == null ? 1000 :
                longestFeatureMap.containsKey(chr) ? longestFeatureMap.get(chr) : 1;
    }

    public void setLongestFeatureMap(Map<String, Integer> longestFeatureMap) {
        this.longestFeatureMap = longestFeatureMap;
    }
}
