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


package org.broad.igv.data.rnai;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.Globals;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;

import java.util.*;

/**
 * @author jrobinso
 */
public class RNAIDataSource implements DataSource {

    /**
     * The screen name (or batch id)
     */
    private String screen;

    /**
     * The screen "condition".  This value can be null.
     */
    private String condition;

    private String displayName;

    private boolean scoresAreSorted = false;

    private Genome genome;

    /**
     * Map of chr -> sorted list of data points.  Data is sorted by increasing
     * start location.
     */
    Map<String, List<LocusScore>> dataMap;

    /**
     * Constructs ...
     *
     * @param screen
     * @param condition
     */
    public RNAIDataSource(String screen, String condition, Genome genome) {
        this.genome = genome;
        this.screen = screen;
        this.condition = condition;
        this.displayName = screen;
        if (condition != null && condition.length() > 0) {
            displayName += " (" + condition + " )";
        }
        dataMap = new HashMap();
    }

    /**
     * Method description
     *
     * @param dpt
     */
    public void addGeneScore(RNAIGeneScore dpt) {
        String chr = dpt.getGene().getChr();
        List<LocusScore> chrDataPoints = dataMap.get(chr);
        if (chrDataPoints == null) {
            chrDataPoints = new ArrayList(500);
            dataMap.put(chr, chrDataPoints);
        }
        chrDataPoints.add(dpt);

        // Also add the score to the "All" chromosome.
        List<LocusScore> chrAllScores = dataMap.get(Globals.CHR_ALL);
        if (chrAllScores == null) {
            chrAllScores = new ArrayList(500);
            dataMap.put(Globals.CHR_ALL, chrAllScores);
        }

        // If a genome is supplied update the "whole genome view"
        if (genome != null) {
            RNAIGeneScore genomeScore = new RNAIGeneScore(dpt);
            int genomeStart = genome.getGenomeCoordinate(dpt.getGene().getChr(), dpt.getStart());
            int genomeEnd = genome.getGenomeCoordinate(dpt.getGene().getChr(), dpt.getEnd());
            genomeScore.setStart(genomeStart);
            genomeScore.setEnd(genomeEnd);
            chrAllScores.add(genomeScore);
        }
    }

    /**
     * Sort gene scores by ascending start position.  This should be called after all gene scores
     * have been added.
     */
    private void sortScores() {
        for (List<LocusScore> scores : dataMap.values()) {
            FeatureUtils.sortFeatureList(scores);
        }

    }


    public double getDataMax() {
        return 3;
    }


    public double getDataMin() {
        return -3;
    }


    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom) {

        if (!scoresAreSorted) {
            sortScores();
        }
        return dataMap.get(chr);
    }

    /**
     * Method description
     *
     * @return
     */
    public TrackType getTrackType() {
        return TrackType.RNAI;
    }

    /**
     * Method description
     *
     * @param statType
     */
    public void setWindowFunction(WindowFunction statType) {

        // Ignored
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isLogNormalized() {
        return true;
    }

    /**
     * Method description
     *
     * @param timestamp
     */
    public void refreshData(long timestamp) {

        // Ignored
    }

    /**
     * The screen name (or batch id)
     *
     * @return
     */
    public String getScreen() {
        return screen;
    }


    /**
     * The screen "condition".  This value can be null.
     *
     * @return
     */
    public String getCondition() {
        return condition;
    }

    public String getName() {
        return displayName;
    }

    /**
     * RNAi data is not windowed
     */
    public WindowFunction getWindowFunction() {
        return null;
    }

    static List<WindowFunction> emptyList = new ArrayList();

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return emptyList;
    }
}
