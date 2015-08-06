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
 * FeatureTrackH5.java
 *
 * Created on November 12, 2007, 8:22 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.*;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SessionXmlAdapters;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.awt.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

/**
 * Represents a track of numeric data
 *
 * @author jrobinso
 */
@XmlType(factoryMethod = "getNextTrack")
public abstract class DataTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(DataTrack.class);

    private DataRenderer renderer;

    // TODO -- memory leak.  This needs to get cleared when the gene list changes
    private HashMap<String, LoadedDataInterval> loadedIntervalCache = new HashMap(200);
    private boolean featuresLoading = false;


    public DataTrack(ResourceLocator locator, String id, String name) {
        super(locator, id, name);
        autoScale = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.CHART_AUTOSCALE);

    }


    public void render(RenderContext context, Rectangle rect){
        if (featuresLoading) {
            return;
        }

        if(isRepeatY(rect)){
            overlay(context, rect);
        }else{
            renderFirstTimeY(context, rect);
            this.lastRenderY = rect.y;
        }
    }


    /**
     * Called the first time we render for a given Y-coordinate
     * @param context
     * @param rect
     */
    private void renderFirstTimeY(RenderContext context, Rectangle rect) {
        List<LocusScore> inViewScores = getInViewScores(context, rect);

        if ((inViewScores == null || inViewScores.size() == 0) && Globals.CHR_ALL.equals(context.getChr())) {
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Data not available for whole genome view; zoom in to see data", rect, g);
        }else{
            getRenderer().render(inViewScores, context, rect, this);
            if(FrameManager.isExomeMode()){
                int x = context.getGraphics().getClipBounds().x;
                Rectangle scaleRect = new Rectangle(x, rect.y, rect.width, rect.height);
                DataRenderer.drawScale(getDataRange(), context, scaleRect);
            }
        }

    }

    public void overlay(RenderContext context, Rectangle rect){
        List<LocusScore> inViewScores = getInViewScores(context, rect);
        if(inViewScores != null){
            synchronized (inViewScores){
                getRenderer().renderScores(this, inViewScores, context, rect);
            }
        }
        getRenderer().renderBorder(this, context, rect);
    }


    private List<LocusScore> getInViewScores(RenderContext context, Rectangle rect){
        String chr = context.getChr();
        int start = (int) context.getOrigin();
        int end = (int) context.getEndLocation() + 1;
        int zoom = context.getZoom();

        List<LocusScore> inViewScores = null;

        LoadedDataInterval interval = loadedIntervalCache.get(context.getReferenceFrame().getName());
        if (interval != null && interval.contains(chr, start, end, zoom)) {
            inViewScores = interval.getScores();
        } else {
            inViewScores = loadScores(context);
        }


        //Not all data sources support whole genome views, tell user if CHR_ALL not available
        if ((inViewScores == null || inViewScores.size() == 0) && Globals.CHR_ALL.equals(chr)) {
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Data not available for whole genome view; zoom in to see data", rect, g);
        } else {
            if (autoScale && !FrameManager.isGeneListMode()) {

                InViewInterval inter = computeScale(start, end, inViewScores);
                if (inter.endIdx > inter.startIdx) {
                    inViewScores = inViewScores.subList(inter.startIdx, inter.endIdx);

                    DataRange dr = getDataRange();
                    float min = Math.min(0, inter.dataMin);
                    float base = Math.max(min, dr.getBaseline());
                    float max = inter.dataMax;
                    // Pathological case where min ~= max  (no data in view)
                    if (max - min <= (2 * Float.MIN_VALUE)) {
                        max = min + 1;
                    }

                    DataRange newDR = new DataRange(min, base, max, dr.isDrawBaseline());
                    newDR.setType(dr.getType());
                    setDataRange(newDR);
                }

            }
        }
        return inViewScores;
    }

    @Override
    public synchronized void load(RenderContext context) {

        String chr = context.getChr();
        int start = (int) context.getOrigin();
        int end = (int) context.getEndLocation() + 1;
        int zoom = context.getZoom();
        LoadedDataInterval interval = loadedIntervalCache.get(context.getReferenceFrame().getName());
        if (interval == null || !interval.contains(chr, start, end, zoom)) {
            loadScores(context);
        }
    }

    public List<LocusScore> loadScores(final RenderContext context) {

        String chr = context.getChr();
        int start = (int) context.getOrigin();
        int end = (int) context.getEndLocation() + 1;
        int zoom = context.getZoom();

        try {
            featuresLoading = true;
            int maxEnd = end;
            Genome genome = GenomeManager.getInstance().getCurrentGenome();

            String queryChr = chr;
            if (genome != null) {
                queryChr = genome.getChromosomeAlias(chr);
                Chromosome c = genome.getChromosome(chr);
                if (c != null) maxEnd = Math.max(c.getLength(), end);
            }

            // Expand interval +/- 50%, unless in a multi-locus mode with "lots" of frames
            boolean multiLocus = FrameManager.isExomeMode() || (FrameManager.getFrames().size() > 4);
            int delta = multiLocus ? 1 : (end - start) / 2;
            int expandedStart = Math.max(0, start - delta);
            int expandedEnd = Math.min(maxEnd, end + delta);
            List<LocusScore> inViewScores = getSummaryScores(queryChr, expandedStart, expandedEnd, zoom);
            LoadedDataInterval interval = new LoadedDataInterval(chr, start, end, zoom, inViewScores);
            loadedIntervalCache.put(context.getReferenceFrame().getName(), interval);
            return inViewScores;

        } finally {
            featuresLoading = false;
        }

    }


    public void clearCaches() {
        loadedIntervalCache.clear();
    }

    public void setRendererClass(Class rc) {
        try {
            renderer = (DataRenderer) rc.newInstance();
        } catch (Exception ex) {
            log.error("Error instantiating renderer ", ex);
        }
    }

    @Override
    protected void setRenderer(Renderer renderer) {
        this.renderer = (DataRenderer) renderer;
    }


    @XmlJavaTypeAdapter(SessionXmlAdapters.Renderer.class)
    @XmlAttribute(name = "renderer")
    @Override
    public DataRenderer getRenderer() {
        if (renderer == null) {
            setRendererClass(getDefaultRendererClass());
        }
        return renderer;
    }


    /**
     * Return a value string for the tooltip window at the given location, or null to signal there is no value
     * at that location
     *
     * @param chr
     * @param position
     * @param y
     * @param frame
     * @return
     */
    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        StringBuffer buf = new StringBuffer();
        LocusScore score = getLocusScoreAt(chr, position, frame);
        // If there is no value here, return null to signal no popup
        if (score == null) {
            return null;
        }
        buf.append(getName() + "<br>");
        if ((getDataRange() != null) && (getRenderer() instanceof XYPlotRenderer)) {
            buf.append("Data scale: " + getDataRange().getMinimum() + " - " + getDataRange().getMaximum() + "<br>");
        }

        buf.append(score.getValueString(position, getWindowFunction()));
        return buf.toString();
    }


    private LocusScore getLocusScoreAt(String chr, double position, ReferenceFrame frame) {
        int zoom = Math.max(0, frame.getZoom());
        List<LocusScore> scores = getSummaryScores(chr, (int) position - 10, (int) position + 10, zoom);


        if (scores == null) {
            return null;
        } else {
            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            int buffer = (int) (2 * bpPerPixel);    /* * */
            return (LocusScore) FeatureUtils.getFeatureAt(position, buffer, scores);
        }
    }


    abstract public List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation, int zoom);

    private InViewInterval computeScale(double origin, double end, List<LocusScore> scores) {

        InViewInterval interval = new InViewInterval();

        if (scores.size() == 1) {
            interval.dataMax = Math.max(0, scores.get(0).getScore());
            interval.dataMin = Math.min(0, scores.get(0).getScore());
        } else {
            interval.startIdx = 0;
            interval.endIdx = scores.size();
            for (int i = 1; i < scores.size(); i++) {
                if (scores.get(i).getEnd() >= origin) {
                    interval.startIdx = i - 1;
                    break;
                }
            }

            for (int i = interval.startIdx + 1; i < scores.size(); i++) {
                LocusScore locusScore = scores.get(i);
                float value = locusScore.getScore();
                if (Float.isNaN(value)) value = 0;
                interval.dataMax = Math.max(interval.dataMax, value);
                interval.dataMin = Math.min(interval.dataMin, value);
                if (locusScore.getStart() > end) {
                    interval.endIdx = i;
                    break;
                }
            }
        }

        return interval;
    }

    /**
     * Get the score over the provided region for the given type. Different types
     * are processed differently. Results are cached according to the provided frameName,
     * if provided. If not, a string is created based on the inputs.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @param type
     * @param frameName
     * @param tracks    Mutation scores require other tracks to calculate the score. If provided,
     *                  use these tracks. If null and not headless we use the currently loaded tracks.
     * @return
     */
    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName,
                                List<Track> tracks) {

        if (end <= start) {
            return 0;
        }
        if (isRegionScoreType(type)) {

            List<LocusScore> scores = null;
            if (frameName == null) {
                //Essentially covering headless case here
                frameName = (chr + start) + end;
                frameName += zoom;
                frameName += type;
            }

            LoadedDataInterval loadedInterval = loadedIntervalCache.get(frameName);
            if (loadedInterval != null && loadedInterval.contains(chr, start, end, zoom)) {
                scores = loadedInterval.getScores();
            } else {
                scores = this.getSummaryScores(chr, start, end, zoom);
                loadedIntervalCache.put(frameName, new LoadedDataInterval(chr, start, end, zoom, scores));
            }


            if (type == RegionScoreType.FLUX) {
                float sumDiffs = 0;
                float lastScore = Float.NaN;
                for (LocusScore score : scores) {
                    if ((score.getEnd() >= start) && (score.getStart() <= end)) {
                        if (Float.isNaN(lastScore)) {
                            lastScore = Math.min(2, Math.max(-2, logScaleData(score.getScore())));
                        } else {
                            float s = Math.min(2, Math.max(-2, logScaleData(score.getScore())));
                            sumDiffs += Math.abs(s - lastScore);
                            lastScore = s;
                        }
                    }
                }
                return sumDiffs;

            } else if (type == RegionScoreType.MUTATION_COUNT) {
                // Sort by overlaid mutation count.
                if (!Globals.isHeadless() && tracks == null) {
                    tracks = IGV.getInstance().getOverlayTracks(this);
                }
                float count = 0;
                String tSamp = this.getSample();
                if (tracks != null && tSamp != null) {
                    for (Track t : tracks) {
                        if (t.getTrackType() == TrackType.MUTATION && tSamp.equals(t.getSample())) {
                            count += t.getRegionScore(chr, start, end, zoom, type, frameName);
                        }
                    }
                }
                return count;
            } else {
                float regionScore = 0;
                int intervalSum = 0;
                boolean hasNan = false;
                for (LocusScore score : scores) {
                    if ((score.getEnd() >= start) && (score.getStart() <= end)) {
                        int interval = Math.min(end, score.getEnd()) - Math.max(start, score.getStart());
                        float value = score.getScore();
                        //For sorting it makes sense to skip NaNs. Not sure about other contexts
                        if(Float.isNaN(value)){
                            hasNan = true;
                            continue;
                        }
                        regionScore += value * interval;
                        intervalSum += interval;
                    }
                }
                if (intervalSum <= 0) {
                    if(hasNan){
                        //If the only existing scores are NaN, the overall score should be NaN
                        return Float.NaN;
                    }else{
                        // No scores in interval
                        return -Float.MAX_VALUE;
                    }
                } else {
                    regionScore /= intervalSum;
                    return (type == RegionScoreType.DELETION) ? -regionScore : regionScore;
                }
            }

        } else {
            return -Float.MAX_VALUE;
        }
    }


    /**
     * Return the average zcore over the interval.
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @return
     */
    public double getAverageScore(String chr, int start, int end, int zoom) {
        double regionScore = 0;
        int intervalSum = 0;
        Collection<LocusScore> scores = getSummaryScores(chr, start, end, zoom);
        for (LocusScore score : scores) {
            if ((score.getEnd() >= start) && (score.getStart() <= end)) {
                int interval = 1; //Math.min(end, score.getEnd()) - Math.max(start, score.getStart());
                float value = score.getScore();
                regionScore += value * interval;
                intervalSum += interval;
            }
        }
        if (intervalSum > 0) {
            regionScore /= intervalSum;
        }
        return regionScore;
    }

    class InViewInterval {
        int startIdx;
        int endIdx;
        float dataMax = 0;
        float dataMin = 0;
    }

    @SubtlyImportant
    private static DataTrack getNextTrack() {
        return (DataTrack) IGVSessionReader.getNextTrack();
    }

}
