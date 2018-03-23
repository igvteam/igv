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

    @Override
    public void dispose() {

    }
}
