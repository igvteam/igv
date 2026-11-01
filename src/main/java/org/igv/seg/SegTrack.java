package org.igv.seg;

import org.igv.feature.FeatureUtils;
import org.igv.feature.Locus;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.feature.mut.MutationRenderer;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.*;
import org.igv.renderer.Renderer;
import org.igv.sample.SampleGroup;
import org.igv.sample.SampleMenuUtils;
import org.igv.sample.SampleSort;
import org.igv.session.RendererFactory;
import org.igv.track.*;
import org.igv.ui.FontManager;
import org.igv.ui.UIConstants;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;
import org.json.JSONObject;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.List;

/**
 * A track with one row per sample.  The type names the content: segmented copy number (seg), mutations (mut),
 * gene expression (gct), or a generic numeric locus-by-sample matrix (igv: .igv, multi-track .tdf, legacy .cn).
 */
public class SegTrack extends AbstractTrack implements ScalableTrack {

    private static final Logger log = LogManager.getLogger(SegTrack.class);

    private static final int EXPANDED_SAMPLE_HEIGHT = 15;
    private static final int SQUISHED_SAMPLE_HEIGHT = 2;

    SegmentedDataSource dataset;
    Renderer renderer;
    TrackType type;
    private transient Rectangle lastClipBounds;
    private final LoadedIntervalCache<Map<String, List<LocusScore>>> loadedIntervalCache = new LoadedIntervalCache<>();

    public SegTrack(ResourceLocator locator,
                    TrackType type,
                    String id,
                    String name,
                    SegmentedDataSource dataset,
                    Genome genome) {

        super(locator, id, name);
        this.defaultExpandedRowHeight = EXPANDED_SAMPLE_HEIGHT;
        this.defaultSquishedRowHeight = SQUISHED_SAMPLE_HEIGHT;
        this.rowHeight = EXPANDED_SAMPLE_HEIGHT;
        this.dataset = dataset;
        this.type = type;
        setDataRange(new DataRange(0, 1, 2));
        setDataType(dataset.getType());
        renderer = switch (type) {
            case mut -> new MutationRenderer();
            case igv, gct -> getDefaultRenderer();
            default -> new HeatmapRenderer();
        };
        initScale(dataset);

        if (type == TrackType.mut) {
            visibilityWindow = 0; // Disable whole-genome view for mutation tracks
        }

        this.initSamples(dataset.getSampleNames());
    }

    public TrackType getType() {
        return type;
    }

    void initScale(SegmentedDataSource dataset) {

        var min = (float) dataset.getDataMin();
        var max = (float) dataset.getDataMax();
        var baseline = 0;

        // If the range is all positive numbers set the min to zero
        if (min > 0) {
            min = 0;
        }

        setDataRange(new DataRange(min, baseline, max));
    }


    @Override
    public int getContentHeight() {
        var nSamples = sampleCount();
        return getSampleHeight() * nSamples + (getSampleGroups().size() - 1) * groupGap;
    }

    @Override
    public int getSampleHeight() {
        return Math.max(1, rowHeight);
    }

    @Override
    public int getNumRows() {
        return sampleCount();
    }

    @Override
    public boolean hasRows() {
        return true;
    }

    @Override
    public void minimizeHeight() {
        setCustomRowHeight(1);
        int newHeight = Math.max(getContentHeight(), getMinimumHeight());
        setHeight(Math.min(newHeight, getHeight()));
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return loadedIntervalCache.contains(frame);
    }

    /**
     * Load the scores of every sample for the frame in one call, expanding the interval as DataTrack does.
     */
    @Override
    public synchronized void load(ReferenceFrame referenceFrame) {

        if (isReadyToPaint(referenceFrame)) return; // already loaded

        Locus expanded = LoadedIntervalCache.expandedInterval(referenceFrame);
        int zoom = referenceFrame.getZoom();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        String queryChr = genome == null ? expanded.chr : genome.getCanonicalChrName(expanded.chr);
        Map<String, List<LocusScore>> features = dataset.getFeatures(queryChr, expanded.start, expanded.end, zoom);
        loadedIntervalCache.put(referenceFrame.getName(), new LoadedDataInterval<>(expanded.chr, expanded.start, expanded.end, zoom, features));
    }

    /**
     * The feature under the position, allowing a tolerance (in bp) on either side so narrow features are not missed.
     * Of the candidates within tolerance, the one nearest the position wins; a feature containing the position has
     * distance 0.  The list must be sorted by start.
     */
    private static LocusScore getFeatureAt(double position, int tolerance, List<LocusScore> scores) {
        if (scores == null || scores.isEmpty()) {
            return null;
        }
        LocusScore closest = null;
        double minDist = Double.POSITIVE_INFINITY;
        for (int i = Math.max(0, FeatureUtils.getIndexBefore(position - tolerance, scores)); i < scores.size(); i++) {
            LocusScore score = scores.get(i);
            if (score.getStart() > position + tolerance) {
                break;
            }
            if (score.getEnd() + tolerance >= position) {
                double dist = Math.max(0, Math.max(score.getStart() - position, position - score.getEnd()));
                if (dist < minDist) {
                    minDist = dist;
                    closest = score;
                }
            }
        }
        return closest;
    }

    /**
     * The loaded scores for the frame, keyed by sample, or null if nothing is loaded for its chromosome.
     */
    private Map<String, List<LocusScore>> getLoadedFeatures(ReferenceFrame frame) {
        LoadedDataInterval<Map<String, List<LocusScore>>> interval = loadedIntervalCache.get(frame.getName());
        return (interval == null || !frame.getChrName().equals(interval.range.chr)) ? null : interval.getFeatures();
    }

    @Override
    public Range getInViewRange(ReferenceFrame frame) {
        Map<String, List<LocusScore>> loaded = getLoadedFeatures(frame);
        if (loaded == null) {
            return null;
        }
        int start = (int) frame.getOrigin();
        int end = (int) frame.getEnd() + 1;
        float min = Float.MAX_VALUE;
        float max = -Float.MAX_VALUE;
        for (List<LocusScore> scores : loaded.values()) {
            if (scores == null) continue;
            for (int i = Math.max(0, FeatureUtils.getIndexBefore(start, scores)); i < scores.size(); i++) {
                LocusScore score = scores.get(i);
                if (score.getStart() > end) break;
                if (score.getEnd() >= start && !Float.isNaN(score.getScore())) {
                    min = Math.min(score.getScore(), min);
                    max = Math.max(score.getScore(), max);
                }
            }
        }
        return min <= max ? new Range(min, max) : null;
    }

    @Override
    public void render(RenderContext context) {

        Rectangle clipBounds = new Rectangle(context.getClipBounds());
        Map<String, List<LocusScore>> loaded = getLoadedFeatures(context.getReferenceFrame());

        final boolean hasGroups = getSampleGroups().size() > 0;
        if (hasGroups && lastClipBounds != null) {
            // Expand the clip bounds to be sure we clear previous labels, but not outside the bounds
            // of the visible rectangle
            Rectangle visibleRect = context.getVisibleRect();
            int expandedTop = Math.max(visibleRect.y, Math.min(clipBounds.y, lastClipBounds.y));
            int expandedBottom = Math.min(visibleRect.y + visibleRect.height,
                    Math.max(clipBounds.y + clipBounds.height, lastClipBounds.y + clipBounds.height));
            clipBounds.y = expandedTop;
            clipBounds.height = expandedBottom - expandedTop;
            context.getGraphics().setClip(clipBounds);
        }

        Rectangle trackRect = context.getTrackRectangle();
        var y = 0;
        for (SampleGroup group : getSampleGroups()) {
            var yLabel = y;

            for (String sample : group.samples()) {
                if (y > clipBounds.y + clipBounds.height) {
                    break;
                }
                if (y + getSampleHeight() > clipBounds.y) {
                    Rectangle r = new Rectangle(trackRect.x, y, trackRect.width, getSampleHeight());
                    var features = loaded == null ? null : loaded.get(sample);
                    this.renderer.render(features, context, r, this);
                }
                y += getSampleHeight();
            }

            String label = group.label();
            if (hasGroups && label != null) {
                var r = new Rectangle(trackRect.x, yLabel, trackRect.width, Math.min(20, y - yLabel));
                context.getGraphics().setColor(UIConstants.getTrackPanelForeground());
                GraphicUtils.drawVerticallyCenteredText(label, 10, r, context.getGraphics(), false, true);
                drawGroupDivider(context.getGraphics(), trackRect, y);
            }

            y += groupGap; // Gap between groups
        }

        lastClipBounds = context.getClipBounds();
    }

    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRect) {

        // Calculate fontsize
        var gap = Math.min(4, getSampleHeight() / 3);
        var fs = Math.min(fontSize, getSampleHeight() - gap);
        var font = FontManager.getFont(fs);
        g2D.setFont(font);
        Rectangle clipRect = g2D.getClipBounds();

        boolean hasGroups = getSampleGroups().size() > 1;
        var y = 0;
        for (var group : getSampleGroups()) {
            for (String s : group.samples()) {
                if (y > clipRect.y + clipRect.height) {
                    break;
                }
                if (y + getSampleHeight() > clipRect.y) {
                    var r = new Rectangle(0, y, visibleRect.width, getSampleHeight());
                    GraphicUtils.drawWrappedText(s, r, g2D, false);
                }
                y += getSampleHeight();
            }
            if (hasGroups) {
                drawGroupDivider(g2D, trackRectangle, y);
            }
            y += groupGap; // Gap between groups
        }
    }

    @Override
    public Renderer<LocusScore> getRenderer() {
        return renderer;
    }

    @Override
    public void setRenderer(Renderer renderer) {
        this.renderer = renderer;
    }

    @Override
    public void setRendererClass(Class rc) {
        try {
            if (Renderer.class.isAssignableFrom(rc)) {
                renderer = (Renderer<LocusScore>) rc.getDeclaredConstructor().newInstance();
            }
        } catch (Exception ex) {
            log.error("Error instantiating renderer ", ex);
        }
    }

    @Override
    public void setWindowFunction(WindowFunction windowFunction) {
        dataset.setWindowFunction(windowFunction);
        loadedIntervalCache.clear();
    }

    @Override
    public WindowFunction getWindowFunction() {
        return dataset.getWindowFunction();
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return dataset.getAvailableWindowFunctions();
    }

    /**
     * The log / non-log (centered on 2) coloring rule is specific to segmented copy number; other types keep the
     * default.
     */
    @Override
    public boolean isLogNormalized() {
        return type != TrackType.seg || dataset.isLogNormalized();
    }

    /**
     * Numeric matrix content, as opposed to segments or mutations.
     */
    @Override
    public boolean isNumeric() {
        return type == TrackType.igv || type == TrackType.gct;
    }

    @Override
    public void unload() {
        super.unload();
        loadedIntervalCache.unsubscribe();
        dataset.dispose();
    }

    @Override
    public List<Component> getPopupMenuItems(final TrackClickEvent te) {

        List<Component> items = new ArrayList<>();

        items.addAll(TrackMenuUtils.getSquishExpandItems(Collections.singletonList(this)));
        items.add(TrackMenuUtils.getRowHeightItem(Collections.singletonList(this)));
        if (isNumeric()) {
            items.add(new JPopupMenu.Separator());
            items.addAll(TrackMenuUtils.getDataMenuItems(Collections.singletonList(this)));
        }

        List<String> keys = AttributeManager.getInstance().getAttributeNames();

        items.add(new JPopupMenu.Separator());
        if (keys.size() > 0) {
            items.add(SampleMenuUtils.getSortByAttributeItem(this));
            items.add(SampleMenuUtils.getGroupByAttributeItem(this));
            items.add(SampleMenuUtils.getFilterByAttributeItem(this));
        }
        items.add(SampleMenuUtils.getFilterByIdItem(this));

        return items;
    }

    /**
     * Return a value string for the tooltip window at the given location, or null to signal there is no value
     * at that location
     *
     * @param chr      The chromosome
     * @param position The genomic position
     * @param mouseX   The mouse X coordinate
     * @param mouseY   Mouse position relative to the enclosing panel
     * @param frame    The reference frame
     * @return A value string for the tooltip, or null.
     */
    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        var trackY = mouseY;
        var y = 0;
        for (var group : getSampleGroups()) {
            var groupPixelHeight = group.samples().size() * getSampleHeight();
            if (trackY >= y && trackY < y + groupPixelHeight) {
                var sampleIndex = (trackY - y) / getSampleHeight();
                var sampleName = group.samples().get(sampleIndex);
                var buf = new StringBuilder();
                var loaded = getLoadedFeatures(frame);
                var score = loaded == null ? null : getFeatureAt(position, (int) (3 * frame.getScale()), loaded.get(sampleName));
                // If there is no value here, return null to signal no popup
                if (score == null) {
                    return null;
                }
                buf.append(sampleName).append("<br>");
                if ((getDataRange() != null) && (getRenderer() instanceof XYPlotRenderer)) {
                    buf.append("Data scale: ").append(getDataRange().getMinimum()).append(" - ").append(getDataRange().getMaximum()).append("<br>");
                }

                buf.append(score.getValueString(position, mouseX, getWindowFunction()));
                return buf.toString();
            } else if (trackY < y + groupPixelHeight + groupGap) {
                return null; // In a gap
            }
            y += groupPixelHeight + groupGap;
        }
        return null;
    }


    /**
     * Get the score over the provided region for the given type. Different types
     * are processed differently. Results are cached according to the provided frameName,
     * if provided. If not, a string is created based on the inputs.
     *
     * @param chr   The chromosome
     * @param start The start of the genomic range
     * @param end   The end of the genomic range
     * @param type  The score type
     */
    @Override
    public void sortSamplesByValue(String chr, int start, int end, RegionScoreType type) {
        sortSamplesByValue(chr, start, end, type == RegionScoreType.DELETION);
    }

    /**
     * Sort samples by their average value over the region, and remember the sort for sessions.  The track's sample
     * list is sorted, as for other sorts, so the order survives regrouping and filtering.  Samples with no value in
     * the region, then samples whose only values are NaN, are sorted to the end in either direction.
     */
    public void sortSamplesByValue(String chr, int start, int end, boolean ascending) {

        if (end <= start || sampleNames == null) {
            return;
        }

        // Compute a value for each sample.  No entry means no scores in the region.  Scores are fetched once for
        // all samples, at the zoom of the current view.
        var sampleScores = new HashMap<String, Float>();
        var allScores = dataset.getFeatures(chr, start, end, FrameManager.getDefaultFrame().getZoom());
        for (String s : sampleNames) {
            var scores = allScores.get(s);
            if (scores == null) {
                continue;
            }
            var regionScore = 0f;
            var intervalSum = 0;
            var hasNan = false;
            for (var score : scores) {
                if ((score.getEnd() >= start) && (score.getStart() <= end)) {
                    var interval = Math.min(end, score.getEnd()) - Math.max(start, score.getStart());
                    var value = score.getScore();
                    //For sorting it makes sense to skip NaNs. Not sure about other contexts
                    if (Float.isNaN(value)) {
                        hasNan = true;
                        continue;
                    }
                    regionScore += value * interval;
                    intervalSum += interval;
                }
            }
            if (intervalSum > 0) {
                sampleScores.put(s, regionScore / intervalSum);
            } else if (hasNan) {
                //If the only existing scores are NaN, the overall score should be NaN
                sampleScores.put(s, Float.NaN);
            }
        }

        sortSamples((o1, o2) -> {
            var v1 = sampleScores.get(o1);
            var v2 = sampleScores.get(o2);
            int rank1 = scoreRank(v1);
            int rank2 = scoreRank(v2);
            if (rank1 != rank2 || rank1 != 0) {
                return Integer.compare(rank1, rank2);
            }
            return ascending ? Float.compare(v1, v2) : Float.compare(v2, v1);
        });

        this.sampleSort = SampleSort.locus(SampleSort.VALUE, chr, start, end, ascending);
    }

    // Values sort first, then no value, then NaN
    private static int scoreRank(Float score) {
        return score == null ? 1 : score.isNaN() ? 2 : 0;
    }

    @Override
    protected void applySampleSort(SampleSort sort) {
        if (sort.getOption() == null || SampleSort.VALUE.equals(sort.getOption())) {    // VALUE is the igv.js default
            sortSamplesByValue(sort.getChr(), sort.getStart(), sort.getEnd(), sort.isAscending());
        } else {
            super.applySampleSort(sort);
        }
    }

    /**
     * The track itself is not filterable, but samples are.
     *
     * @return false
     */
    @Override
    public boolean isFilterable() {
        return false;
    }

    public static boolean isSegFile(String pathOrURL) {
        String path;
        try {
            var url = new URL(pathOrURL);
            path = url.getPath();
        } catch (MalformedURLException e) {
            // Not a valid URL, assume it's a file path
            path = pathOrURL;
        }
        var lowerPath = path.toLowerCase();
        return lowerPath.endsWith(".seg") || lowerPath.endsWith(".seg.gz") || lowerPath.endsWith(".seg.txt") ||
                lowerPath.endsWith(".seg.txt.gz");
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
            json.put("renderer", renderer.getClass().getName());
        }
    }

    @Override
    public void unmarshalJSON(JSONObject jsonObject) {
        super.unmarshalJSON(jsonObject);
        if (jsonObject.has("renderer")) {
            Class rendererClass = RendererFactory.getRendererClass(jsonObject.getString("renderer"));
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
