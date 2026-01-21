/*
* GisticScoreList.java
*
* Created on June 21, 2007, 8:29 AM
*
* To change this template, choose Tools | Template Manager
* and open the template in the editor.
*/
package org.igv.track;


import org.igv.logging.*;
import org.igv.feature.IGVFeature;
import org.igv.feature.FeatureUtils;
import org.igv.feature.GisticScore;
import org.igv.feature.LocusScore;
import org.igv.renderer.DataRange;
import org.igv.renderer.GisticTrackRenderer;
import org.igv.renderer.Renderer;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public class GisticTrack extends AbstractTrack {

    private static Logger log = LogManager.getLogger(GisticTrack.class);

    private static final int DEFAULT_HEIGHT = 50;

    private double maxQValue = 0;

    private double maxGScore = 0;

    private boolean hasScores = false;

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

    public GisticTrack() {
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return hasScores;
    }

    @Override
    public void load(ReferenceFrame referenceFrame) {
        //
    }

    @Override
    public int getMinimumHeight() {
        return 25;
    }


    public void setScores(List<GisticScore> scores) {

        this.hasScores = true;

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
     * @param mouseX
     *@param ignore  @return
     */
    public String getValueStringAt(String chr, double position, int mouseX, int ignore, ReferenceFrame frame) {

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

    // GisticTrack does not expose its renderer
    public Renderer getRenderer() {
        return null;
    }

}
