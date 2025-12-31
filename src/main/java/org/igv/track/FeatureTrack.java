package org.igv.track;

import htsjdk.tribble.Feature;
import htsjdk.tribble.NamedFeature;
import htsjdk.tribble.TribbleException;
import org.igv.event.IGVEvent;
import org.igv.logging.*;
import org.igv.Globals;
import org.igv.event.DataLoadedEvent;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
import org.igv.feature.*;
import org.igv.feature.Range;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.*;
import org.igv.tools.motiffinder.MotifFinderSource;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.MessageUtils;
import org.igv.util.BrowserLauncher;
import org.igv.util.ResourceLocator;
import org.igv.util.StringUtils;
import org.igv.variant.VariantTrack;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

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
    protected static final Color SELECTED_FEATURE_ROW_COLOR = new Color(100, 100, 100, 30);
    private static final int DEFAULT_EXPANDED_HEIGHT = 35;
    private static final int DEFAULT_SQUISHED_HEIGHT = 12;

    private int expandedRowHeight = DEFAULT_EXPANDED_HEIGHT;
    private int squishedRowHeight = DEFAULT_SQUISHED_HEIGHT;

    boolean fatalLoadError = false;

    Track.DisplayMode lastFeatureMode = null;  // Keeps track of the feature display mode before an auto-switch to COLLAPSE

    // TODO -- this is a memory leak, this cache needs cleared when the reference frame collection (gene list) changes
    /**
     * Map of reference frame name -> packed features
     */
    protected Map<String, PackedFeatures<Feature>> packedFeaturesMap = Collections.synchronizedMap(new HashMap<>());

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

    private static boolean drawBorder = true;

    private boolean alternateExonColor = false;

    private String trackLine = null;

    private boolean groupByStrand = false;
    private String labelField;

    public FeatureTrack() {

    }

    // TODO -- there are WAY too many constructors for this class

    /**
     * Constructor used by SpliceJunctionTrack, BlatTrack, and MotifTrack
     *
     * @param locator -- For splice junctions, ResourceLocator for associated alignment track.  Null otherwise
     * @param id
     * @param name
     */
    public FeatureTrack(ResourceLocator locator, String id, String name) {
        super(locator, id, name);
        setSortable(false);
    }

    /**
     * Constructor with no ResourceLocator.  Note:  tracks using this constructor will not be recorded in the
     * "Resources" section of session files.  This method is used by ".genome" and gbk genome loaders.
     */
    public FeatureTrack(String id, String name, FeatureSource source) {
        super(null, id, name);
        init(null, source);
        setSortable(false);
    }

    /**
     * Constructor specifically for BigWig data source
     *
     * @param locator
     * @param id
     * @param name
     * @param source
     */
    public FeatureTrack(ResourceLocator locator, String id, String name, FeatureSource source) {
        super(locator, id, name);
        init(locator, source);
        setSortable(false);
    }

    public FeatureTrack(ResourceLocator locator, FeatureSource source) {
        super(locator);
        init(locator, source);
        setSortable(false);
    }


    /**
     * SMap files (bionano), MutationTrack
     */
    public FeatureTrack(ResourceLocator locator, String id, FeatureSource source) {
        super(locator, id, locator.getTrackName());
        init(locator, source);
    }

    /**
     * Create a new track which is a shallow copy of this one.  Currently used by SashimiPlot.
     */
    public FeatureTrack(FeatureTrack featureTrack) {
        this(featureTrack.getId(), featureTrack.getName(), featureTrack.source);
    }

    protected void init(ResourceLocator locator, FeatureSource source) {

        this.source = source;
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

        IGVEventBus.getInstance().subscribe(DataLoadedEvent.class, this);

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
        if (e instanceof DataLoadedEvent) {
//            DataLoadedEvent event = (DataLoadedEvent) e;
//            if (IGV.hasInstance()) {
//                // TODO -- WHY IS THIS HERE????
//                //TODO Assuming this is necessary, there can be many data loaded events in succession,
//                //don't want to layout for each one
//                IGV.getInstance().layoutMainPanel();
//            }
        } else {
            log.warn("Unknown event type: " + e.getClass());
        }
    }

    @Override
    public int getHeight() {
        if (!isVisible()) {
            return 0;
        }
        int rowHeight = getDisplayMode() == DisplayMode.SQUISHED ? squishedRowHeight : expandedRowHeight;
        int minHeight = margin + rowHeight * Math.max(1, getNumberOfFeatureLevels());
        return Math.max(minHeight, super.getHeight());
    }

    public int getExpandedRowHeight() {
        return expandedRowHeight;
    }

    public void setExpandedRowHeight(int expandedRowHeight) {
        this.expandedRowHeight = expandedRowHeight;
    }

    public int getSquishedRowHeight() {
        return squishedRowHeight;
    }

    public void setSquishedRowHeight(int squishedRowHeight) {
        this.squishedRowHeight = squishedRowHeight;
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
                if (scoreType == RegionScoreType.MUTATION_COUNT && this.getTrackType() == TrackType.MUTATION) {
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
            int maxNumber =  10;
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
                        buf.append("<hr><br<b>+ " + (allFeatures.size() - maxNumber) +  " more</b>");
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

        int rowHeight;
        DisplayMode mode = getDisplayMode();
        switch (mode) {
            case SQUISHED:
                rowHeight = getSquishedRowHeight();
                break;
            case EXPANDED:
                rowHeight = getExpandedRowHeight();
                break;
            default:
                rowHeight = getHeight();
        }

        return Math.max(0, (y - this.getY() - this.margin) / rowHeight);

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

        PackedFeatures<Feature> packedFeatures = packedFeaturesMap.get(frame.getName());

        if (packedFeatures == null) {
            return null;
        }

        List<PackedFeatures<Feature>.FeatureRow> rows = packedFeatures.getRows();
        if (featureRow < 0 || featureRow >= rows.size()) {
            return null;
        }

        //If features are stacked we look at only the row.
        //If they are collapsed on top of each other, we get all features in all rows
        List<Feature> possFeatures = rows.get(featureRow).getFeatures();


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
            if (i == selectedFeatureRowIndex)
                setSelectedFeatureRowIndex(FeatureTrack.NO_FEATURE_ROW_SELECTED);
            else {
                //make this track selected
                setSelected(true);
                //select the appropriate row
                setSelectedFeatureRowIndex(i);
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


    public void overlay(RenderContext context, Rectangle rect) {
        renderFeatures(context, rect);
    }

    @Override
    public void setDisplayMode(DisplayMode mode) {
        // Explicity setting the display mode overrides the automatic switch
        lastFeatureMode = null;

        for (PackedFeatures pf : packedFeaturesMap.values()) {
            pf.pack(mode, groupByStrand);
        }

        super.setDisplayMode(mode);
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        if (!isShowFeatures(frame)) {
            packedFeaturesMap.clear();
            return true;  // Ready by definition (nothing to paint)
        } else {
            PackedFeatures packedFeatures = packedFeaturesMap.get(frame.getName());
            String chr = frame.getChrName();
            int start = (int) frame.getOrigin();
            int end = (int) frame.getEnd();
            return (packedFeatures != null && packedFeatures.containsInterval(chr, start, end));
        }
    }

    public void load(ReferenceFrame frame) {
        loadFeatures(frame.getChrName(), (int) frame.getOrigin(), (int) frame.getEnd(), frame.getName());
    }

    /**
     * Loads and segregates features into rows such that they do not overlap.
     *
     * @param chr
     * @param start
     * @param end
     */
    protected void loadFeatures(final String chr, final int start, final int end, final String frame) {

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
    public void render(RenderContext context, Rectangle rect) {
        Rectangle renderRect = new Rectangle(rect);
        renderRect.y = renderRect.y + margin;
        renderRect.height -= margin;


        showFeatures = isShowFeatures(context.getReferenceFrame());
        if (showFeatures) {
            if (lastFeatureMode != null) {
                super.setDisplayMode(lastFeatureMode);
                lastFeatureMode = null;
            }
            renderFeatures(context, renderRect);
        } else if (coverageRenderer != null) {
            if (getDisplayMode() != DisplayMode.COLLAPSED) {
                // An ugly hack, but we want to prevent this for vcf tracks
                if (!(this instanceof VariantTrack)) {
                    lastFeatureMode = getDisplayMode();
                    super.setDisplayMode(DisplayMode.COLLAPSED);
                }
            }
            renderCoverage(context, renderRect);
        }

        if (FeatureTrack.drawBorder) {
            Graphics2D borderGraphics = context.getGraphic2DForColor(UIConstants.TRACK_BORDER_GRAY);
            borderGraphics.drawLine(rect.x, rect.y, rect.x + rect.width, rect.y);
            //TODO Fix height for variant track
            borderGraphics.drawLine(rect.x, rect.y + rect.height, rect.x + rect.width, rect.y + rect.height);
        }

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

        if (scores == null) {
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            Rectangle textRect = new Rectangle(inputRect);

            // Keep text near the top of the track rectangle
            textRect.height = Math.min(inputRect.height, 20);
            String message = getZoomInMessage(chr);
            GraphicUtils.drawCenteredText(message, textRect, g);

        } else {
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
     * @param inputRect
     */
    protected void renderFeatures(RenderContext context, Rectangle inputRect) {

        if (log.isTraceEnabled()) {
            String msg = String.format("renderFeatures: %s frame: %s", getName(), context.getReferenceFrame().getName());
            log.trace(msg);
        }

        PackedFeatures packedFeatures = packedFeaturesMap.get(context.getReferenceFrame().getName());

        if (packedFeatures == null || !packedFeatures.overlapsInterval(context.getChr(), (int) context.getOrigin(), (int) context.getEndLocation() + 1)) {
            return;
        }

        try {
            renderFeatureImpl(context, inputRect, packedFeatures);
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

    protected void renderFeatureImpl(RenderContext context, Rectangle inputRect, PackedFeatures packedFeatures) {


        Renderer renderer = getRenderer();

        if (getDisplayMode() == DisplayMode.COLLAPSED && !isGroupByStrand()) {
            List<Feature> features = packedFeatures.getFeatures();
            if (features != null) {
                renderer.render(features, context, inputRect, this);
            }
        } else {
            List<PackedFeatures.FeatureRow> rows = packedFeatures.getRows();
            if (rows != null && rows.size() > 0) {

                // Divide rectangle into equal height levels
                double h = getDisplayMode() == DisplayMode.SQUISHED ? squishedRowHeight : expandedRowHeight;
                Rectangle rect = new Rectangle(inputRect.x, inputRect.y, inputRect.width, (int) h);
                int i = 0;

                if (renderer instanceof FeatureRenderer) ((FeatureRenderer) renderer).reset();
                for (PackedFeatures.FeatureRow row : rows) {

                    renderer.render(row.features, context, new Rectangle(rect), this);

                    if (selectedFeatureRowIndex == i) {
                        Graphics2D fontGraphics = context.getGraphic2DForColor(SELECTED_FEATURE_ROW_COLOR);
                        fontGraphics.fillRect(rect.x, rect.y, rect.width, rect.height);
                    }

                    rect.y += h;
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
        PackedFeatures packedFeatures = packedFeaturesMap.get(frame.getName());

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

    public static boolean isDrawBorder() {
        return drawBorder;
    }

    public static void setDrawBorder(boolean drawBorder) {
        FeatureTrack.drawBorder = drawBorder;
    }

    public boolean isAlternateExonColor() {
        return alternateExonColor;
    }

    /**
     * This method exists for Plugin tracks. When restoring a session there is no guarantee of track
     * order, so arguments referring to other tracks may fail to resolve. We update those references
     * here after all tracks have been processed
     *
     * @param allTracks
     */
    public void updateTrackReferences(List<Track> allTracks) {

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
        PackedFeatures<Feature> packedFeatures = packedFeaturesMap.get(frame.getName());
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


    @Override
    public void setTrackLine(String trackLine) {
        this.trackLine = trackLine;
    }

    /**
     * Return "track" line information for exporting features to a file.  Default is null, subclasses may override.
     *
     * @return
     */
    public String getExportTrackLine() {
        return trackLine;
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
    public void marshalXML(Document document, Element element) {
        element.setAttribute("groupByStrand", String.valueOf(groupByStrand));
        if (labelField != null) {
            element.setAttribute("featureNameProperty", labelField);
        }
        super.marshalXML(document, element);

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

}

