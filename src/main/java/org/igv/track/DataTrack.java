/*
 * FeatureTrackH5.java
 *
 * Created on November 12, 2007, 8:22 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.igv.track;

import org.igv.Globals;
import org.igv.feature.FeatureUtils;
import org.igv.feature.LocusScore;
import org.igv.feature.Locus;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.DataRenderer;
import org.igv.renderer.GraphicUtils;
import org.igv.renderer.Renderer;
import org.igv.renderer.XYPlotRenderer;
import org.igv.session.RendererFactory;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;
import org.json.JSONObject;
import org.w3c.dom.Element;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * Represents a track of numeric data
 *
 * @author jrobinso
 */

public abstract class DataTrack extends AbstractTrack implements ScalableTrack {

    public static final Color ZOOM_IN_TEXT_COLOR = Globals.isDarkMode() ? Color.lightGray.brighter() : Color.gray;
    private static Logger log = LogManager.getLogger(DataTrack.class);

    private DataRenderer renderer;

    private LoadedIntervalCache<List<LocusScore>> loadedIntervalCache = new LoadedIntervalCache<>();

    public DataTrack(ResourceLocator locator, String id, String name) {
        super(locator, id, name);
        setHeight(PreferencesManager.getPreferences().getAsInt(Constants.CHART_TRACK_HEIGHT_KEY));
    }

    public DataTrack() {
    }

    @Override
    public TrackType getType() {
        return TrackType.wig;
    }

    @Override
    public List<Component> getPopupMenuItems(TrackClickEvent te) {
        return TrackMenuUtils.getDataMenuItems(Collections.singleton(this));
    }

    @Override
    public int getContentHeight() {
        return height;
    }

    @Override
    public boolean isNumeric() {
        return true;
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return loadedIntervalCache.contains(frame);
    }


    public synchronized void load(ReferenceFrame referenceFrame) {

        if (isReadyToPaint(referenceFrame)) return; // already loaded

        Locus expanded = LoadedIntervalCache.expandedInterval(referenceFrame);
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        String queryChr = genome == null ? expanded.chr : genome.getCanonicalChrName(expanded.chr);
        LoadedDataInterval<List<LocusScore>> interval = getSummaryScores(queryChr, expanded.start, expanded.end, referenceFrame.getZoom());
        loadedIntervalCache.put(referenceFrame.getName(), interval);
    }


    public void render(RenderContext context) {

        List<LocusScore> inViewScores = getInViewScores(context.getReferenceFrame());

        if ((inViewScores == null || inViewScores.size() == 0) && Globals.CHR_ALL.equals(context.getChr())) {
            Graphics2D g = context.getGraphic2DForColor(ZOOM_IN_TEXT_COLOR);
            Rectangle visibleRect = context.getVisibleRect();
            GraphicUtils.drawCenteredText("Data not available for whole genome view; zoom in to see data", visibleRect, g);
        } else {
            getRenderer().render(inViewScores, context, context.getTrackRectangle(), this);
        }
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

    @Override
    public void unload() {
        super.unload();
        loadedIntervalCache.unsubscribe();
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
                float count = 0;
                String tSamp = this.getSample();
                if (tracks != null && tSamp != null) {
                    for (Track t : tracks) {
                        if (t.getDataType() == DataType.MUTATION && tSamp.equals(t.getSample())) {
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

    @Override
    public void marshalJSON(JSONObject json) {

        super.marshalJSON(json);

        if (renderer != null) {
            RendererFactory.RendererType type = RendererFactory.getRenderType(renderer);
            if (type != null) {
                json.put("graphType", type.name());
            }
        }
    }

    @Override
    public void unmarshalJSON(JSONObject jsonObject) {

        super.unmarshalJSON(jsonObject);

        if (jsonObject.has("graphType")) {

            Class rendererClass = RendererFactory.getRendererClass(jsonObject.getString("graphType"));

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
