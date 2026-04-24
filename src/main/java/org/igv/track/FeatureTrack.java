package org.igv.track;

import htsjdk.tribble.Feature;
import htsjdk.tribble.NamedFeature;
import htsjdk.tribble.TribbleException;
import org.igv.Globals;
import org.igv.event.IGVEvent;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
import org.igv.event.ViewChange;
import org.igv.feature.*;
import org.igv.feature.Range;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.*;
import org.igv.renderer.Renderer;
import org.igv.tools.motiffinder.MotifFinderSource;
import org.igv.ui.IGV;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.MessageUtils;
import org.igv.util.BrowserLauncher;
import org.igv.util.ResourceLocator;
import org.igv.util.StringUtils;
import org.igv.variant.VariantTrack;
import org.json.JSONObject;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

import static org.igv.track.TrackMenuUtils.*;

/**
 * Track which displays features, typically showing regions of the genome
 * in a qualitative way. Features are rendered using the specified FeatureRenderer.
 * The gene track is an example of a feature track.
 *
 * @author jrobinso
 */


public class FeatureTrack extends AbstractTrack implements IGVEventObserver {

    private static Logger log = LogManager.getLogger(FeatureTrack.class);


    public static final int MINIMUM_FEATURE_SPACING = 0;
    public static final int DEFAULT_MARGIN = 5;
    public static final int NO_FEATURE_ROW_SELECTED = -1;
    private static final int DEFAULT_EXPANDED_HEIGHT = 35;
    private static final int DEFAULT_SQUISHED_HEIGHT = 12;
    private int squishedRowHeight = DEFAULT_SQUISHED_HEIGHT;
    private transient int maxFeatureRow = 1;

    boolean fatalLoadError = false;
    /**
     * Map of reference frame -> packed features
     */
    protected Map<ReferenceFrame, PackedFeatures<PackedFeature>> packedFeaturesMap = Collections.synchronizedMap(new HashMap<>());

    protected Renderer renderer;

    private DataRenderer coverageRenderer;

    // true == features,  false =  coverage
    private boolean showFeatures = true;

    public FeatureSource source;

    //track which row of the expanded track is selected by the user. Selection goes away if tracks are collpased
    protected int selectedFeatureRowIndex = NO_FEATURE_ROW_SELECTED;

    //Feature selected by the user.  This is repopulated on each handleDataClick() call.
    protected IGVFeature selectedFeature = null;

    int margin = DEFAULT_MARGIN;

    private boolean alternateExonColor = false;

    private boolean groupByStrand = false;
    private String labelField;

    public FeatureTrack() {

    }

    @Override
    public TrackType getType() {
        return TrackType.annotation;
    }

    @Override
    public List<Component> getPopupMenuItems(TrackClickEvent te) {

        Collection<Track> tracks = Collections.singleton(this);

        List<Component> items = new ArrayList<>();

        for (Component item : getDisplayModeMenuItems(tracks)) {
            items.add(item);
        }
        items.add(new JSeparator());

        items.add(getGroupByStrandItem(tracks));

        if (tracks.size() == 1) {
            Track t = tracks.iterator().next();
            Feature f = t.getFeatureAtMousePosition(te);

            ReferenceFrame frame = te.getFrame();
            if (frame == null && !FrameManager.isGeneListMode()) {
                frame = FrameManager.getDefaultFrame();
            }

            String featureName = "";
            if (f != null) {
                items.add(new JPopupMenu.Separator());
                items.add(getCopyDetailsItem(f, te));

                Feature sequenceFeature = f;
                if (sequenceFeature instanceof IGVFeature) {
                    featureName = ((IGVFeature) sequenceFeature).getName();
                    double position = te.getChromosomePosition();
                    Collection<Exon> exons = ((IGVFeature) sequenceFeature).getExons();
                    if (exons != null) {
                        for (Exon exon : exons) {
                            if (position > exon.getStart() && position < exon.getEnd()) {
                                sequenceFeature = exon;
                                break;
                            }
                        }
                    }
                }

                items.add(getCopySequenceItem(sequenceFeature));

                if (frame != null && PreferencesManager.getPreferences().get(Constants.EXTVIEW_URL) != null) {
                    Range r = frame.getCurrentRange();
                    items.add(getExtendViewItem(featureName, sequenceFeature, r));
                }

                items.add(getBlatItem(sequenceFeature));
            }
        }

        items.add(new JPopupMenu.Separator());
        items.add(getChangeFeatureWindow(tracks));

        items.add(new JPopupMenu.Separator());
        items.add(getShowFeatureNames(tracks));
        items.add(getFeatureNameAttribute(tracks));

        return items;

    }

    public FeatureTrack(ResourceLocator locator, FeatureSource source) {
        super(locator);
        init(locator, source);
    }

    /**
     * Constructor used by SpliceJunctionTrack, BlatTrack, and MotifTrack
     *
     * @param locator -- For splice junctions, ResourceLocator for associated alignment track.  Null otherwise
     * @param id
     * @param name
     */
    public FeatureTrack(ResourceLocator locator, String id, String name) {
        super(locator, id, name);
        this.rowHeight = DEFAULT_EXPANDED_HEIGHT;
    }

    /**
     * Constructor with no ResourceLocator.  Used by ".genome" and gbk genome loaders.
     */
    public FeatureTrack(String id, String name, FeatureSource source) {
        super(null, id, name);
        init(null, source);
    }


    /**
     * Create a new track which is a shallow copy of this one.  Currently used by SashimiPlot.
     */
    public FeatureTrack(FeatureTrack featureTrack) {
        super(null, featureTrack.getId(), featureTrack.getName());
        init(null, featureTrack.source);
    }

    protected void init(ResourceLocator locator, FeatureSource source) {

        this.source = source;

        this.rowHeight = DEFAULT_EXPANDED_HEIGHT;
        setMinimumHeight(10);

        coverageRenderer = new BarChartRenderer();

        Integer vizWindow = locator == null ? null : locator.getVisibilityWindow();

        if (vizWindow != null) {
            visibilityWindow = vizWindow;
        } else {
            int sourceFeatureWindowSize = source.getFeatureWindowSize();
            int defVisibilityWindow = PreferencesManager.getPreferences().getAsInt(Constants.DEFAULT_VISIBILITY_WINDOW);
            if (sourceFeatureWindowSize > 0 && defVisibilityWindow > 0) {  // Only apply a default if the feature source supports visibility window.
                visibilityWindow = defVisibilityWindow * 1000;
            } else {
                visibilityWindow = sourceFeatureWindowSize;
            }
        }

        this.renderer = locator != null && locator.getPath().endsWith("junctions.bed") ?
                new SpliceJunctionRenderer() : new IGVFeatureRenderer();

        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(ViewChange.class, this);

    }

    @Override
    public void unload() {
        super.unload();
        if (source != null) {
            source.close();
        }
    }

    /**
     * Called after features are finished loading, which can be asynchronous
     */
    public void receiveEvent(IGVEvent e) {
        if (e instanceof FrameManager.ChangeEvent) {
            // Build a set of current frames
            Set<ReferenceFrame> currentFrames = new HashSet<>(FrameManager.getFrames());

            // Remove cached intervals for frames that are no longer present
            // packedFeaturesMap is a synchronizedMap; synchronize while iterating/removing
            synchronized (packedFeaturesMap) {
                packedFeaturesMap.keySet().removeIf(frame -> !currentFrames.contains(frame));
            }
        } else if (e instanceof ViewChange) {
            if (!((ViewChange) e).panning) {
                updateMaxFeatureRow();
            }
        } else {
            log.warn("Unknown event type: " + e.getClass());
        }
    }

    @Override
    public boolean hasDisplayMode() {
        return true;
    }

    @Override
    public int getContentHeight() {
        if (!isVisible()) {
            return 0;
        }
        int minHeight = margin + Math.round(getRowHeight() * Math.max(1, maxFeatureRow));
        return Math.max(minHeight, super.getContentHeight());
    }

    @Override
    public int getNumRows() {
        return Math.max(1, maxFeatureRow);
    }

    private void updateMaxFeatureRow() {
        maxFeatureRow = 1;
        if (getDisplayMode() != DisplayMode.COLLAPSED) {
            if (packedFeaturesMap.size() > 0) {
                for (ReferenceFrame frame : packedFeaturesMap.keySet()) {
                    PackedFeatures pf = packedFeaturesMap.get(frame);
                    maxFeatureRow = Math.max(maxFeatureRow, pf.getMaxPackedRow((int) frame.getOrigin(), (int) frame.getEnd()));
                }
            }
        }
    }

    public void setRendererClass(Class rc) {
        try {
            renderer = (Renderer) rc.getDeclaredConstructor().newInstance();
        } catch (Exception ex) {
            log.error("Error instatiating renderer ", ex);
        }
    }

    public void setMargin(int margin) {
        this.margin = margin;
    }

    @Override
    public void setProperties(TrackProperties trackProperties) {
        super.setProperties(trackProperties);
        if (trackProperties.getFeatureVisibilityWindow() >= 0) {
            setVisibilityWindow(trackProperties.getFeatureVisibilityWindow());
        }
        alternateExonColor = trackProperties.isAlternateExonColor();

    }


    public void setWindowFunction(WindowFunction type) {
        // Ignored for feature tracks
    }


    /**
     * Return the maximum number of features for any panel in this track.  In whole genome view there is a single panel,
     * but there are multiple in gene list view (one for each gene list).
     *
     * @return
     */
    public int getNumberOfFeatureLevels() {
        if (packedFeaturesMap.size() > 0) {
            int n = 1;
            for (PackedFeatures pf : packedFeaturesMap.values()) {
                //dhmay adding null check.  To my mind this shouldn't be necessary, but we're encountering
                //it intermittently.  Food for future thought
                if (pf != null) {
                    n = Math.max(n, pf.getRowCount());
                }
            }
            return n;
        }
        return 1;
    }

    /**
     * Return a score over the interval.  This is required by the track interface to support sorting.
     */
    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType scoreType, String frameName) {

        try {
            Iterator<Feature> features = source.getFeatures(chr, start, end);
            if (features != null) {
                if (scoreType == RegionScoreType.MUTATION_COUNT && this.getDataType() == DataType.MUTATION) {
                    int count = 0;
                    while (features.hasNext()) {
                        Feature f = features.next();
                        if (f.getStart() > end) {
                            break;
                        }
                        if (f.getEnd() >= start) {
                            count++;
                        }
                    }

                    return count;
                } else if (scoreType == RegionScoreType.SCORE) {
                    // Average score of features in region.  Note: Should the score be weighted by genomic size?
                    float regionScore = 0;
                    int nValues = 0;
                    while (features.hasNext()) {
                        Feature f = features.next();
                        if (f instanceof IGVFeature) {
                            if ((f.getEnd() >= start) && (f.getStart() <= end)) {
                                float value = ((IGVFeature) f).getScore();
                                regionScore += value;
                                nValues++;
                            }
                        }
                    }
                    if (nValues == 0) {
                        // No scores in interval
                        return -Float.MAX_VALUE;
                    } else {
                        return regionScore / nValues;
                    }
                }
            }
        } catch (IOException e) {
            log.error("Error counting features.", e);
        }
        return -Float.MAX_VALUE;
    }


    public Renderer getRenderer() {
        if (renderer == null) {
            setRenderer(new IGVFeatureRenderer());
        }
        return renderer;
    }

    /**
     * Return a string for popup text.
     *
     * @param chr
     * @param position in genomic coordinates
     * @param mouseX
     * @return
     */
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        if (showFeatures) {

            List<Feature> allFeatures = getAllFeaturesContaining(position, mouseY, frame);
            if (allFeatures == null) {
                return null;
            }

            StringBuffer buf = new StringBuffer();
            boolean firstFeature = true;
            int maxNumber = 10;
            int n = 1;
            for (Feature feature : allFeatures) {
                if (feature != null && feature instanceof IGVFeature) {
                    if (!firstFeature) {
                        buf.append("<hr><br>");
                    }

                    IGVFeature igvFeature = (IGVFeature) feature;
                    String vs = igvFeature.getValueString(position, mouseX, null);
                    buf.append(vs);

                    if (IGV.getInstance().isShowDetailsOnClick()) {
                        String url = getFeatureURL(igvFeature);
                        if (url != null) {
                            buf.append("<br/><a href=\"" + url + "\">" + url + "</a>");
                        }
                    }
                    firstFeature = false;
                    if (n > maxNumber) {
                        buf.append("<hr><br<b>+ " + (allFeatures.size() - maxNumber) + " more</b>");
                        break;
                    }
                }
                n++;
            }
            return buf.toString();
        } else {
            int zoom = Math.max(0, frame.getZoom());
            if (source == null) {
                return null;
            }
            List<LocusScore> scores = source.getCoverageScores(chr, (int) position - 10, (int) position + 10, zoom);

            if (scores == null) {
                return "";
            } else {
                // give a +/- 2 pixel buffer, otherwise very narrow features will be missed.
                double bpPerPixel = frame.getScale();
                int minWidth = (int) (2 * bpPerPixel);    /* * */
                LocusScore score = (LocusScore) FeatureUtils.getFeatureAt(position, minWidth, scores);
                return score == null ? null : "Mean count: " + score.getScore();
            }

        }
    }


    /**
     * Return an info URL for a specific feature, constructed as follows
     *
     * @param igvFeature
     * @return
     */
    private String getFeatureURL(IGVFeature igvFeature) {
        String url = igvFeature.getURL();    // Explicity URL setting
        if (url == null) {
            String trackURL = getFeatureInfoURL();   // Template
            if (trackURL != null) {
                String idOrName = igvFeature.getIdentifier() != null ?
                        igvFeature.getIdentifier() :
                        labelField != null ?
                                igvFeature.getDisplayName(labelField) :
                                igvFeature.getName();
                url = trackURL.replaceAll("\\$\\$", StringUtils.encodeURL(idOrName));
            }
        }
        return url;
    }

    @Override
    public String getLabelField() {
        return labelField;
    }

    public void setLabelField(String labelField) {
        this.labelField = labelField;
    }

    /**
     * Get all features which overlap the specified locus
     *
     * @return
     */
    public List<Feature> getFeatures(String chr, int start, int end) {
        List<Feature> features = new ArrayList<Feature>();
        try {
            Iterator<Feature> iter = source.getFeatures(chr, start, end);
            while (iter.hasNext()) {
                Feature f = iter.next();
                if (f.getEnd() >= start && f.getStart() <= end) {
                    features.add(f);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return features;
    }

    /**
     * @param position in genomic coordinates
     * @param y        pixel location in panel coordinates.  // TODO offset by track origin before getting here?
     * @param frame
     * @return
     */
    protected List<Feature> getAllFeaturesContaining(double position, int y, ReferenceFrame frame) {
        // Determine the level number (for expanded tracks)
        int featureRow = getFeatureRow(y);
        return getFeaturesAtPositionInFeatureRow(position, featureRow, frame);
    }

    /**
     * Determine which row the user clicked in and return the appropriate feature
     *
     * @param y
     * @return
     */
    private int getFeatureRow(int y) {
        return Math.max(0, (int) (y / getRowHeight()));
    }

    /**
     * Knowing the feature row, figure out which feature is at {@code position}.
     *
     * @param position
     * @param featureRow
     * @param frame
     * @return
     */
    public List<Feature> getFeaturesAtPositionInFeatureRow(double position, int featureRow, ReferenceFrame frame) {

        PackedFeatures<PackedFeature> packedFeatures = packedFeaturesMap.get(frame);

        if (packedFeatures == null) {
            return null;
        }

        List<PackedFeatures<PackedFeature>.FeatureRow> rows = packedFeatures.getRows();
        if (featureRow < 0 || featureRow >= rows.size()) {
            return null;
        }

        //If features are stacked we look at only the row.
        //If they are collapsed on top of each other, we get all features in all rows
        List<PackedFeature> possFeatures = rows.get(featureRow).getFeatures();


        List<Feature> featureList = null;
        if (possFeatures != null) {
            // give a minum 2 pixel or 1/2 bp window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            double flanking = 4 * bpPerPixel;
            featureList = FeatureUtils.getAllFeaturesContaining(position, flanking, possFeatures);
        }
        return featureList;
    }


    public WindowFunction getWindowFunction() {
        return WindowFunction.count;
    }

    @Override
    public boolean handleDataClick(TrackClickEvent te) {

        MouseEvent e = te.getMouseEvent();

        // Toggle selection of row
        if (getDisplayMode() != DisplayMode.COLLAPSED) {
            int i = getFeatureRow(e.getY());
            if (i == selectedFeatureRowIndex) {
                setSelectedFeatureRowIndex(FeatureTrack.NO_FEATURE_ROW_SELECTED);
            }
        }


        //For feature selection
        selectedFeature = null;

        Feature f = getFeatureAtMousePosition(te);
        if (f != null && f instanceof IGVFeature) {
            IGVFeature igvFeature = (IGVFeature) f;
            if (selectedFeature != null && igvFeature.contains(selectedFeature) && (selectedFeature.contains(igvFeature))) {
                //If something already selected, then if it's the same as this feature, deselect, otherwise, select
                //this feature.
                //todo: contains() might not do everything I want it to.
                selectedFeature = null;
            } else {
                //if nothing already selected, or something else selected,
                // select this feature
                selectedFeature = igvFeature;
            }

            if (IGV.getInstance().isShowDetailsOnClick()) {
                openTooltipWindow(te);
            } else {
                String url = getFeatureURL(igvFeature);
                if (url != null) {
                    try {
                        BrowserLauncher.openURL(url);
                    } catch (IOException e1) {
                        log.error("Error launching url: " + url);
                    }
                    e.consume();
                    return true;
                }
            }
        }

        return false;
    }

    public Feature getFeatureAtMousePosition(TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();
        final ReferenceFrame referenceFrame = te.getFrame();
        if (referenceFrame != null) {
            double location = referenceFrame.getChromosomePosition(e);
            List<Feature> features = getAllFeaturesContaining(location, e.getY(), referenceFrame);
            return (features != null && features.size() > 0) ? features.get(0) : null;
        } else {
            return null;
        }
    }

    /**
     * Required by the interface, really not applicable to feature tracks
     */
    public boolean isLogNormalized() {
        return true;
    }


    @Override
    public void setDisplayMode(DisplayMode mode) {

        // Deal with the legacy "squished" mode.  This is an expanded mode with a reduced row height
        if (mode == DisplayMode.SQUISHED) {
            setRowHeight(DEFAULT_SQUISHED_HEIGHT);
            mode = DisplayMode.EXPANDED;
        }

        super.setDisplayMode(mode);

        for (PackedFeatures pf : packedFeaturesMap.values()) {
            pf.pack(mode, groupByStrand);
        }
        updateMaxFeatureRow();

    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        if (!isShowFeatures(frame)) {
            packedFeaturesMap.clear();
            return true;  // Ready by definition (nothing to paint)
        } else {
            PackedFeatures packedFeatures = packedFeaturesMap.get(frame);
            String chr = frame.getChrName();
            int start = (int) frame.getOrigin();
            int end = (int) frame.getEnd();
            return (packedFeatures != null && packedFeatures.containsInterval(chr, start, end));
        }
    }

    public void load(ReferenceFrame frame) {
        loadFeatures(frame.getChrName(), (int) frame.getOrigin(), (int) frame.getEnd(), frame);
        updateMaxFeatureRow();
    }

    /**
     * Loads and segregates features into rows such that they do not overlap.
     *
     * @param chr
     * @param start
     * @param end
     * @param frame
     */
    protected void loadFeatures(final String chr, final int start, final int end, final ReferenceFrame frame) {

        if (source == null) {
            return;
        }

        try {
            int expandedStart;
            int expandedEnd;
            int vw = getVisibilityWindow();
            if (vw > 0) {
                int delta = (end - start) / 2;
                expandedStart = start - delta;
                expandedEnd = end + delta;
                if (expandedEnd < 0) {
                    expandedEnd = Integer.MAX_VALUE;  // overflow
                }
            } else {
                expandedStart = 0;
                expandedEnd = Integer.MAX_VALUE;
            }


            //Make sure we are only querying within the chromosome we allow for somewhat pathological cases of start
            //being negative and end being outside, but only if directly queried. Our expansion should not
            //set start < 0 or end > chromosomeLength
            if (start >= 0) {
                expandedStart = Math.max(0, expandedStart);
            }

            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            if (genome != null) {
                Chromosome c = genome.getChromosome(chr);
                if (c != null && end < c.getLength()) expandedEnd = Math.min(c.getLength(), expandedEnd);
            }

            Iterator<Feature> iter = source.getFeatures(chr, expandedStart, expandedEnd);

            if (iter == null) {
                PackedFeatures pf = new PackedFeatures(chr, expandedStart, expandedEnd);
                packedFeaturesMap.put(frame, pf);
            } else {
                PackedFeatures pf = new PackedFeatures(chr, expandedStart, expandedEnd, iter, this.getDisplayMode(), groupByStrand);
                packedFeaturesMap.put(frame, pf);
                //log.warn("Loaded " + chr + " " + expandedStart + "-" + expandedEnd);
            }

        } catch (Exception e) {
            // Mark the interval with an empty feature list to prevent an endless loop of load attempts.
            PackedFeatures pf = new PackedFeatures(chr, start, end);
            packedFeaturesMap.put(frame, pf);
            String msg = "Error loading features for interval: " + chr + ":" + start + "-" + end + " <br>" + e.toString();
            MessageUtils.showMessage(msg);
            log.error(msg, e);
        }

    }

    @Override
    public void render(RenderContext context) {

        // Draw entire track.
        Rectangle renderRect = context.getTrackRectangle();

        context.getReferenceFrame().getCurrentRange();

        renderRect.y = renderRect.y + margin;
        renderRect.height -= margin;

        showFeatures = isShowFeatures(context.getReferenceFrame());
        if (showFeatures) {
            renderFeatures(context, renderRect);
        } else if (coverageRenderer != null) {
            renderCoverage(context, renderRect);
        }
    }

    @Override
    public boolean supportsWholeGenome() {
        return source == null || source.supportsWholeGenome();
    }

    protected boolean isShowFeatures(ReferenceFrame frame) {

        if (frame.getChrName().equals(Globals.CHR_ALL)) {
            return false;
        } else {
            double windowSize = frame.getEnd() - frame.getOrigin();
            int vw = getVisibilityWindow();
            return (vw <= 0 || windowSize <= vw);
        }
    }

    protected void renderCoverage(RenderContext context, Rectangle inputRect) {

        final String chr = context.getChr();
        List<LocusScore> scores = (source != null && chr.equals(Globals.CHR_ALL)) ?
                source.getCoverageScores(chr, (int) context.getOrigin(),
                        (int) context.getEndLocation(), context.getZoom()) :
                null;

        if (scores != null) {
            float max = getMaxEstimate(scores);
            if (this.dataRange == null) {
                setDataRange(new DataRange(0, 0, max));
            } else {
                this.dataRange.maximum = max;
            }
            coverageRenderer.render(scores, context, inputRect, this);
        }
    }

    protected String getZoomInMessage(String chr) {
        return chr.equals(Globals.CHR_ALL) ? "Select a chromosome and zoom in to see features." :
                "Zoom in to see features, or right-click to increase Feature Visibility Window.";
    }

    private float getMaxEstimate(List<LocusScore> scores) {
        float max = 0;
        int n = scores.size();
        for (int i = 0; i < n; i++) {
            max = Math.max(max, scores.get(i).getScore());
        }
        return max;
    }

    /**
     * Render features in the given input rectangle.
     *
     * @param context
     * @param ignore  - deprecated
     */
    protected void renderFeatures(RenderContext context, Rectangle ignore) {

        if (log.isTraceEnabled()) {
            String msg = String.format("renderFeatures: %s frame: %s", getName(), context.getReferenceFrame().getName());
            log.trace(msg);
        }

        PackedFeatures packedFeatures = packedFeaturesMap.get(context.getReferenceFrame());

        if (packedFeatures == null || !packedFeatures.overlapsInterval(context.getChr(), (int) context.getOrigin(), (int) context.getEndLocation() + 1)) {
            return;
        }

        try {
            renderFeatureImpl(context, packedFeatures);
        } catch (TribbleException e) {
            log.error("Tribble error", e);
            //Error loading features.  We'll let the user decide if this is "fatal" or not.
            if (!fatalLoadError) {
                fatalLoadError = true;
                boolean unload = MessageUtils.confirm("<html> Error loading features: " + e.getMessage() +
                        "<br>Unload track " + getName() + "?");
                if (unload) {
                    Collection<Track> tmp = Arrays.asList((Track) this);
                    IGV.getInstance().deleteTracks(tmp);
                    IGV.getInstance().repaint();
                } else {
                    fatalLoadError = false;
                }
            }
        }
    }

    protected void renderFeatureImpl(RenderContext context, PackedFeatures packedFeatures) {

        Rectangle trackRectangle = context.getTrackRectangle();
        Rectangle clipBounds = context.getClipBounds();

        Renderer renderer = getRenderer();

        if (getDisplayMode() == DisplayMode.COLLAPSED && !isGroupByStrand()) {
            List<Feature> features = packedFeatures.getFeatures();
            if (features != null) {
                renderer.render(features, context, trackRectangle, this);
            }
        } else {
            List<PackedFeatures.FeatureRow> rows = packedFeatures.getRows();
            if (rows != null && rows.size() > 0) {

                // Divide rectangle into equal height levels
                int intH = Math.max(1, Math.round(getRowHeight()));
                Rectangle rect = new Rectangle(trackRectangle.x, trackRectangle.y, trackRectangle.width, intH);
                int i = 0;

                if (renderer instanceof FeatureRenderer) ((FeatureRenderer) renderer).reset();
                for (PackedFeatures.FeatureRow row : rows) {

                    if (rect.y > clipBounds.y + clipBounds.height) {
                        break;  // Gone past clip bounds
                    }
                    if (rect.y + rect.height < clipBounds.y) {
                        // Not yet in clip bounds
                        rect.y += intH;
                        i++;
                        continue;
                    }

                    renderer.render(row.features, context, rect, this);

                    //  if (selectedFeatureRowIndex == i) {
                    //     Graphics2D fontGraphics = context.getGraphic2DForColor(SELECTED_FEATURE_ROW_COLOR);
                    //      fontGraphics.fillRect(rect.x, rect.y, rect.width, rect.height);
                    //  }

                    rect.y += intH;
                    i++;
                }
            }
        }
    }


    /**
     * Return the nextLine or previous feature relative to the center location.
     * TODO -- investigate delegating this method to FeatureSource, where it might be possible to simplify the implementation
     *
     * @param chr
     * @param center
     * @param forward
     * @return
     * @throws IOException
     */
    public Feature nextFeature(String chr, double center, boolean forward, ReferenceFrame frame) throws IOException {

        Feature f = null;
        boolean canScroll = (forward && !frame.windowAtEnd()) || (!forward && frame.getOrigin() > 0);
        PackedFeatures packedFeatures = packedFeaturesMap.get(frame);

        // Compute a buffer to define "next"
        double buffer = Math.max(1, frame.getScale());

        if (packedFeatures != null && packedFeatures.containsInterval(chr, (int) center - 1, (int) center + 1)) {
            if (packedFeatures.getFeatures().size() > 0 && canScroll) {
                List<Feature> centerSortedFeatures = packedFeatures.getCenterSortedFeatures();
                f = (forward ?
                        FeatureUtils.getFeatureCenteredAfter(center + buffer, centerSortedFeatures) :
                        FeatureUtils.getFeatureCenteredBefore(center - buffer, centerSortedFeatures));
            }
        }
        if (f == null) {
            int searchBuferSize = (int) (frame.getScale() * 1000);
            f = FeatureTrackUtils.nextFeature(source, chr, packedFeatures.getStart(), packedFeatures.getEnd(), center, searchBuferSize, forward);
        }


        return f;
    }


    public void setVisibilityWindow(int windowSize) {
        super.setVisibilityWindow(windowSize);
        packedFeaturesMap.clear();
    }

    public int getSelectedFeatureRowIndex() {
        return selectedFeatureRowIndex;
    }

    public void setSelectedFeatureRowIndex(int selectedFeatureRowIndex) {
        this.selectedFeatureRowIndex = selectedFeatureRowIndex;
    }

    public IGVFeature getSelectedFeature() {
        return selectedFeature;
    }

    public boolean isAlternateExonColor() {
        return alternateExonColor;
    }

    /**
     * Features are packed upon loading, effectively a cache.
     * This clears that cache. Used to force a refresh
     *
     * @api
     */
    public void clearPackedFeatures() {
        this.packedFeaturesMap.clear();
    }

    /**
     * Return currently loaded features.  Used to export features to a file.
     *
     * @param frame
     * @return
     */
    public List<Feature> getVisibleFeatures(ReferenceFrame frame) {
        PackedFeatures<PackedFeature> packedFeatures = packedFeaturesMap.get(frame);
        if (packedFeatures == null) {
            return Collections.emptyList();
        } else {
            Range currentRange = frame.getCurrentRange();
            return packedFeatures.getFeatures()
                    .stream().filter(igvFeature ->
                            igvFeature.getEnd() >= currentRange.getStart() &&
                                    igvFeature.getStart() <= currentRange.getEnd() &&
                                    igvFeature.getChr().equals(currentRange.getChr()))
                    .collect(Collectors.toList());
        }
    }


    public void setGroupByStrand(boolean selected) {
        this.groupByStrand = selected;
        for (PackedFeatures pf : packedFeaturesMap.values()) {
            pf.pack(getDisplayMode(), groupByStrand);
        }
    }

    public boolean isGroupByStrand() {
        return groupByStrand;
    }

    @Override
    public boolean isSearchable() {
        return source != null && source.isSearchable();
    }

    @Override
    public List<NamedFeature> search(String token) {
        return source.search(token);
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("featureNameProperty")) {
            this.labelField = element.getAttribute("featureNameProperty");
        }

        this.groupByStrand = "true".equals(element.getAttribute("groupByStrand"));

        NodeList tmp = element.getElementsByTagName("SequenceMatchSource");
        if (tmp.getLength() > 0) {
            Element sourceElement = (Element) tmp.item(0);
            MotifFinderSource source = new MotifFinderSource();
            source.unmarshalXML(sourceElement, version);
            this.source = source;
        }
    }

    @Override
    public void marshalJSON(JSONObject jsonObject) {
        super.marshalJSON(jsonObject);
        if (labelField != null) {
            jsonObject.put("nameField", labelField);
        }
        if (groupByStrand) {
            jsonObject.put("groupBy", "strand");
        }
    }

    @Override
    public void unmarshalJSON(JSONObject jsonObject) {

        super.unmarshalJSON(jsonObject);

        if (jsonObject.has("nameField")) {
            this.labelField = jsonObject.getString("nameField");
        }
        if (jsonObject.has("groupBy")) {
            this.groupByStrand = "strand".equals(jsonObject.getString("groupBy"));
        }
    }
}

