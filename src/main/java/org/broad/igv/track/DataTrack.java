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

import org.broad.igv.event.IGVEvent;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.DataRenderer;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.renderer.XYPlotRenderer;
import org.broad.igv.session.RendererFactory;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * Represents a track of numeric data
 *
 * @author jrobinso
 */

public abstract class DataTrack extends AbstractTrack implements ScalableTrack, IGVEventObserver {

    private static Logger log = LogManager.getLogger(DataTrack.class);

    private DataRenderer renderer;

    private Map<String, LoadedDataInterval<List<LocusScore>>> loadedIntervalCache = new HashMap(200);

    public DataTrack(ResourceLocator locator, String id, String name) {
        super(locator, id, name);
        loadedIntervalCache = Collections.synchronizedMap(new HashMap<>());
        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
    }

    public DataTrack() {
    }

    public void receiveEvent(IGVEvent event) {

        if (event instanceof FrameManager.ChangeEvent) {

            Collection<ReferenceFrame> frames = ((FrameManager.ChangeEvent) event).frames();
            Map<String, LoadedDataInterval<List<LocusScore>>> newCache = Collections.synchronizedMap(new HashMap<>());
            for (ReferenceFrame f : frames) {
                newCache.put(f.getName(), loadedIntervalCache.get(f.getName()));
            }
            loadedIntervalCache = newCache;


        } else {
            log.warn("Unknown event type: " + event.getClass());
        }
    }

    @Override
    public boolean isNumeric() {
        return true;
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        String chr = frame.getChrName();
        int start = (int) frame.getOrigin();
        int end = (int) frame.getEnd();
        int zoom = frame.getZoom();
        LoadedDataInterval interval = loadedIntervalCache.get(frame.getName());
        return (interval != null && interval.contains(chr, start, end, zoom));

    }


    public synchronized void load(ReferenceFrame referenceFrame) {

        if (isReadyToPaint(referenceFrame)) return; // already loaded

        String chr = referenceFrame.getChrName();
        int start = (int) referenceFrame.getOrigin();
        int end = (int) referenceFrame.getEnd() + 1;
        int zoom = referenceFrame.getZoom();
        int maxEnd = end;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        String queryChr = chr;
        if (genome != null) {
            queryChr = genome.getCanonicalChrName(chr);
            Chromosome c = genome.getChromosome(chr);
            if (c != null) maxEnd = Math.max(c.getLength(), end);
        }

        // Expand interval +/- 50%, unless in a multi-locus mode with "lots" of frames
        boolean multiLocus = (FrameManager.getFrames().size() > 4);
        int delta = multiLocus ? 1 : (end - start) / 2;
        int expandedStart = Math.max(0, start - delta);
        int expandedEnd = Math.min(maxEnd, end + delta);
        if(expandedEnd < 0) {
            // overflow
            expandedEnd = Integer.MAX_VALUE;
        }
        LoadedDataInterval<List<LocusScore>> interval = getSummaryScores(queryChr, expandedStart, expandedEnd, zoom);
        loadedIntervalCache.put(referenceFrame.getName(), interval);

    }


    public void render(RenderContext context, Rectangle rect) {

        List<LocusScore> inViewScores = getInViewScores(context.getReferenceFrame());

        if ((inViewScores == null || inViewScores.size() == 0) && Globals.CHR_ALL.equals(context.getChr())) {
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Data not available for whole genome view; zoom in to see data", rect, g);
        } else {
            getRenderer().render(inViewScores, context, rect, this);
        }

    }


    public void overlay(RenderContext context, Rectangle rect) {

        List<LocusScore> inViewScores = getInViewScores(context.getReferenceFrame());

        if ((inViewScores == null || inViewScores.size() == 0) && Globals.CHR_ALL.equals(context.getChr())) {
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Data not available for whole genome view; select chromosome to see data", rect, g);
        } else if (inViewScores != null) {
            synchronized (inViewScores) {
                getRenderer().renderScores(this, inViewScores, context, rect);
            }
        }
        getRenderer().renderBorder(this, context, rect);
    }


    public List<LocusScore> getInViewScores(ReferenceFrame referenceFrame) {

        LoadedDataInterval<List<LocusScore>> interval = loadedIntervalCache.get(referenceFrame.getName());
        String chr = referenceFrame.getChrName();
        int start = (int) referenceFrame.getOrigin();
        int end = (int) referenceFrame.getEnd() + 1;
        int zoom = referenceFrame.getZoom();
        if (interval == null || !(chr.equals(interval.range.chr))) { // Try the data we have, even if not perfect !interval.contains(chr, start, end, zoom)) {
            return Collections.EMPTY_LIST;
        }

        List<LocusScore> inViewScores = interval.getFeatures();

        // Trim scores
        int startIdx = Math.max(0, FeatureUtils.getIndexBefore(start, inViewScores));
        int endIdx = inViewScores.size();   // Starting guess
        int tmp = FeatureUtils.getIndexBefore(end, inViewScores);

        if (tmp < 0)
            return Collections.EMPTY_LIST;

        else {
            for (int i = tmp; i < inViewScores.size(); i++) {
                if (inViewScores.get(i).getStart() > end) {
                    endIdx = i + 1;
                    break;
                }
            }
            endIdx = Math.max(startIdx + 1, endIdx);

            return startIdx == 0 && endIdx == inViewScores.size() ?
                    inViewScores :
                    inViewScores.subList(startIdx, endIdx);
        }
    }

    public Range getInViewRange(ReferenceFrame referenceFrame) {

        List<LocusScore> scores = getInViewScores(referenceFrame);
        if (scores.size() > 0) {
            float min = Float.MAX_VALUE;
            float max = -Float.MAX_VALUE;
            for (LocusScore score : scores) {
                float value = score.getScore();
                if (!Float.isNaN(value)) {
                    min = Math.min(value, min);
                    max = Math.max(value, max);
                }
            }
            return new Range(min, max);
        } else {
            return null;
        }

    }


    public void clearCaches() {
        loadedIntervalCache.clear();
    }

    public void setRendererClass(Class rc) {
        try {
            renderer = (DataRenderer) rc.getDeclaredConstructor().newInstance();
        } catch (Exception ex) {
            log.error("Error instantiating renderer ", ex);
        }
    }

    @Override
    public void setRenderer(Renderer renderer) {
        this.renderer = (DataRenderer) renderer;
    }


    @Override
    public DataRenderer getRenderer() {
        if (renderer == null) {
            renderer = (DataRenderer) getDefaultRenderer();
        }
        return renderer;
    }


    /**
     * Return a value string for the tooltip window at the given location, or null to signal there is no value
     * at that location
     *
     * @param chr
     * @param position
     * @param mouseX
     * @param frame    @return
     */
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {
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

        buf.append(score.getValueString(position, mouseX, getWindowFunction()));
        return buf.toString();
    }


    private LocusScore getLocusScoreAt(String chr, double position, ReferenceFrame frame) {
        int zoom = Math.max(0, frame.getZoom());
        List<LocusScore> scores = getSummaryScores(chr, (int) position - 10, (int) position + 10, zoom).getFeatures();


        if (scores == null) {
            return null;
        } else {
            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            int buffer = (int) (2 * bpPerPixel);    /* * */
            return (LocusScore) FeatureUtils.getFeatureAt(position, buffer, scores);
        }
    }


    abstract public LoadedDataInterval<List<LocusScore>> getSummaryScores(String chr, int startLocation, int endLocation, int zoom);

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

            LoadedDataInterval<List<LocusScore>> loadedInterval = loadedIntervalCache.get(frameName);
            if (loadedInterval != null && loadedInterval.contains(chr, start, end, zoom)) {
                scores = loadedInterval.getFeatures();
            } else {
                scores = this.getSummaryScores(chr, start, end, zoom).getFeatures();
                loadedIntervalCache.put(frameName, new LoadedDataInterval<>(chr, start, end, zoom, scores));
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
                        if (Float.isNaN(value)) {
                            hasNan = true;
                            continue;
                        }
                        regionScore += value * interval;
                        intervalSum += interval;
                    }
                }
                if (intervalSum <= 0) {
                    if (hasNan) {
                        //If the only existing scores are NaN, the overall score should be NaN
                        return Float.NaN;
                    } else {
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
        Collection<LocusScore> scores = getSummaryScores(chr, start, end, zoom).getFeatures();
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

    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        if (renderer != null) {
            RendererFactory.RendererType type = RendererFactory.getRenderType(renderer);
            if (type != null) {
                element.setAttribute("renderer", type.name());
            }
        }

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("renderer")) {

            Class rendererClass = RendererFactory.getRendererClass(element.getAttribute("renderer"));
            if (rendererClass != null) {
                try {
                    renderer = (DataRenderer) rendererClass.newInstance();

                } catch (Exception e) {
                    log.error("Error instantiating renderer: " + rendererClass.getName(), e);
                }
            }
        }

    }

}
