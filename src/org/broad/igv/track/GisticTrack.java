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
* GisticScoreList.java
*
* Created on June 21, 2007, 8:29 AM
*
* To change this template, choose Tools | Template Manager
* and open the template in the editor.
*/
package org.broad.igv.track;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.GisticScore;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.GisticTrackRenderer;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class GisticTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(GisticTrack.class);

    private static final int DEFAULT_HEIGHT = 50;

    private double maxQValue = 0;

    private double maxGScore = 0;

    /**
     * Map of chromosome -> sorted list of scores
     */
    Map<String, List<GisticScore>> ampScoreMap;

    Map<String, List<GisticScore>> delScoreMap;

    GisticTrackRenderer renderer;


    public GisticTrack(ResourceLocator locator) {
        super(locator);
        ampScoreMap = new HashMap<String, List<GisticScore>>();
        delScoreMap = new HashMap<String, List<GisticScore>>();
        renderer = new GisticTrackRenderer();
        setHeight(DEFAULT_HEIGHT);
        renderer = new GisticTrackRenderer();
        setSortable(false);

    }

    @Override
    public int getMinimumHeight() {
        return 25;
    }


    /**
     * Method description
     *
     * @param scores
     */
    public void setScores(List<GisticScore> scores) {

        for (GisticScore score : scores) {
            String chr = score.getChromosome();
            if (score.getType() == GisticScore.Type.AMP) {
                addScore(score, chr, ampScoreMap);
            } else {
                addScore(score, chr, delScoreMap);

            }
        }
        updateMaxValues(scores);
    }

    protected void addScore(GisticScore score, String chr, Map<String, List<GisticScore>> map) {
        List<GisticScore> scoreList = map.get(chr);
        if (scoreList == null) {
            scoreList = new ArrayList();
            map.put(chr, scoreList);
        }
        scoreList.add(score);
    }

    private void updateMaxValues(List<GisticScore> scores) {
        for (GisticScore score : scores) {
            if (!Double.isInfinite(maxGScore) && score.getGScore() > maxGScore) {
                maxGScore = score.getGScore();
            }
            if (!Double.isInfinite(maxQValue) && score.getQValue() > maxQValue) {
                maxQValue = score.getQValue();
            }
        }
        setDataRange(new DataRange(0, 0, (float) maxQValue));
    }

    /**
     * Method description
     *
     * @param chr
     * @return
     */
    public List<GisticScore> getAmpScores(String chr) {
        return ampScoreMap.get(chr);
    }

    /**
     * Method description
     *
     * @param chr
     * @return
     */
    public List<GisticScore> getDelScores(String chr) {
        return delScoreMap.get(chr);
    }

    /**
     * Return the score with the maximum "G Score" over the specified region.
     * Assumes the scores are sorted by location
     *
     * @param chr
     * @param start
     * @param end
     * @return
     */
    public GisticScore getMaxForRegion(String chr, long start, long end) {

        List<GisticScore> scores = ampScoreMap.get(chr);
        if ((scores == null) || scores.isEmpty()) {
            return null;
        }

        GisticScore maxScore = null;
        for (GisticScore score : scores) {
            if (maxScore == null) {
                if ((score.getStart() >= start)
                        || ((start >= score.getStart()) && (start <= score.getEnd()))) {
                    maxScore = score;
                }
            } else if (score.getStart() > end) {
                break;
            } else {
                if (score.getGScore() > maxScore.getGScore()) {
                    maxScore = score;
                }
            }
        }

        // If we haven't found bounding scores yet the region is to the right off the last
        // score.  Use the last score for the region.  This should be a rare case and should
        // probably be logged.
        return (maxScore == null) ? ((GisticScore) scores.get(scores.size() - 1)) : maxScore;

    }

    /**
     * Method description
     *
     * @return
     */
    public double getMaxQValue() {
        return maxQValue;
    }

    /**
     * Method description
     *
     * @return
     */
    public double getMaxGScore() {
        return maxGScore;
    }

    /**
     * Method description
     *
     * @param context
     * @param rect
     */
    public void render(RenderContext context, Rectangle rect) {
        if (renderer == null) {
            log.error("Null renderer !!");

        } else {
            renderer.render(this, context, rect);
        }
    }

    double dataMax = -1;

    /**
     * Method description
     *
     * @return
     */
    public double getDataMax() {
        if (dataMax < 0) {
            dataMax = getMaxGScore();
        }
        return dataMax;
    }

    /**
     * Method description
     *
     * @param type
     */
    public void setWindowFunction(WindowFunction type) {

        // ignored
    }

    /**
     * Method description
     *
     * @return
     */
    public WindowFunction getWindowFunction() {
        return WindowFunction.median;
    }

    /**
     * Method description
     *
     * @param chr
     * @return
     */
    public double getMedian(String chr) {
        return 1.0;
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
     * This method is required for the interface, but will not be called.
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @return
     */
    public List<IGVFeature> getFeatures(String chr, int startLocation, int endLocation) {
        return null;
    }

    /**
     * This method is required for the interface, but will not be called.
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @param zoom
     * @return
     */
    public List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation,
                                             int zoom) {

        List<LocusScore> ss = new ArrayList();
        List<GisticScore> ampScores = ampScoreMap.get(chr);
        if (ampScores != null) {
            ss.addAll(ampScores);
        }
        List<GisticScore> delScores = delScoreMap.get(chr);
        if (delScores != null) {
            ss.addAll(delScores);
        }
        return ss;

        // return getFeatures(chr, startLocation, endLocation);
    }

    /**
     * This method is required for the interface, but will not be called.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @param type
     * @param frameName
     * @return
     */
    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0;
    }


    /**
     * Method description
     *
     * @param chr
     * @param position
     * @param ignore
     * @return
     */
    public String getValueStringAt(String chr, double position, int ignore, ReferenceFrame frame) {

        // give a 2 pixel window, otherwise very narrow features will be missed.
        double bpPerPixel = frame.getScale();
        double minWidth = 2 * bpPerPixel;    /*
                                              */

        LocusScore amp = null;
        List<GisticScore> ampScores = ampScoreMap.get(chr);
        if (ampScores != null) {
            amp = (LocusScore) FeatureUtils.getFeatureAt(position, 0, ampScores);
        }

        LocusScore del = null;
        List<GisticScore> delScores = delScoreMap.get(chr);
        if (delScores != null) {
            del = (LocusScore) FeatureUtils.getFeatureAt(position, 0, delScores);
        }

        if ((amp == null) && (del == null)) {
            return "";
        } else {
            StringBuffer buf = new StringBuffer();
            if (amp != null) {
                buf.append("Amplification score: " + amp.getScore());
            }
            if (del != null) {
                buf.append("<br>Deletion score: " + del.getScore());
            }
            return buf.toString();

        }
    }

    /**
     * Method description
     *
     * @param chr
     * @param startLocation
     * @param endLocation
     * @return
     */
    public List<List<IGVFeature>> getFeaturesByLevels(String chr, int startLocation, int endLocation) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Method description
     *
     * @return
     */
    public int getNumberOfFeatureLevels() {
        return 1;
    }

    /**
     * Method description
     *
     * @param isMultiLevel
     */
    public void setMultiLevelFeatures(boolean isMultiLevel) {

        // ignore
    }

    /**
     * Method description
     *
     * @return
     */
    public boolean isMuliLevelFeatures() {
        throw new UnsupportedOperationException("Not supported yet.");
    }


    public Color getColor() {
        return null;
    }

    public void setColor(Color color) {
        // Not used
    }

    public Color getAltColor() {
        return null;
    }

    public void setAltColor(Color color) {
        // Not used
    }

    // GisticTrack does not expose its renderer
    public Renderer getRenderer() {
        return null;
    }


}
