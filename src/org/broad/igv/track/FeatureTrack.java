/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */
package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.igv.util.collections.CachedIntervals;
import org.broad.igv.variant.VariantTrack;
import org.broad.tribble.Feature;
import org.broad.tribble.TribbleException;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class FeatureTrack extends AbstractTrack {

    private static Logger log = Logger.getLogger(FeatureTrack.class);

    public static final int MINIMUM_FEATURE_SPACING = 5;
    public static final int DEFAULT_MARGIN = 5;
    public static final int NO_FEATURE_ROW_SELECTED = -1;
    protected static final Color SELECTED_FEATURE_ROW_COLOR = new Color(100, 100, 100, 30);
    private static final int DEFAULT_EXPANDED_HEIGHT = 35;
    private static final int DEFAULT_SQUISHED_HEIGHT = 12;

    private int expandedRowHeight = DEFAULT_EXPANDED_HEIGHT;
    private int squishedRowHeight = DEFAULT_SQUISHED_HEIGHT;

    boolean fatalLoadError = false;

    Track.DisplayMode lastFeatureMode = null;  // Keeps track of the feature display mode before an auto-switch to COLLAPSE


    protected List<Rectangle> levelRects = new ArrayList();

    // TODO -- this is a memory leak, this cache needs cleared when the reference frame collection (gene list) changes
    /**
     * Map of reference frame name -> packed features
     */
    protected Map<String, PackedFeatures<IGVFeature>> packedFeaturesMap = new HashMap();

    private FeatureRenderer renderer = new IGVFeatureRenderer();

    private DataRenderer coverageRenderer;

    // true == features,  false =  coverage
    private boolean showFeatures = true;

    protected FeatureSource source;

    protected boolean featuresLoading = false;

    //track which row of the expanded track is selected by the user.
    //Selection goes away if tracks are collpased
    protected int selectedFeatureRowIndex = NO_FEATURE_ROW_SELECTED;

    //Feature selected by the user.  This is repopulated on each handleDataClick() call.
    protected IGVFeature selectedFeature = null;


    int margin = DEFAULT_MARGIN;

    private static boolean drawBorder = true;

    private boolean alternateExonColor = false;

    /**
     * Construct with no feature source.  Currently this is only used for the SpliceJunctionFinderTrack subclass.
     *
     * @param id
     * @param name
     */
    public FeatureTrack(String id, String name) {
        super(id, name);
        setSortable(false);
    }

    public FeatureTrack(ResourceLocator locator, String id, String name) {
        super(locator, id, name);
        setSortable(false);
    }

    /**
     * Constructor with no ResourceLocator.  Note:  tracks using this constructor will not be recorded in the
     * "Resources" section of session files.   Currently only used for the default gene track.
     *
     * @param id
     * @param name
     * @param source
     */
    public FeatureTrack(String id, String name, FeatureSource source) {
        super(id, name);
        init(source);
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
        init(source);
        setSortable(false);
    }


    public FeatureTrack(ResourceLocator locator, FeatureSource source) {
        super(locator);
        init(source);
        this.getMinimumHeight();
        setSortable(false);
    }


    public FeatureTrack(ResourceLocator locator, String id, FeatureSource source) {
        super(locator, id);
        init(source);
    }

    protected void init(FeatureSource source) {

        this.source = source;
        setMinimumHeight(10);
        setColor(Color.blue.darker());

        coverageRenderer = new BarChartRenderer();
        if (source.getFeatureWindowSize() > 0) {
            visibilityWindow = source.getFeatureWindowSize();
        }
    }

    @Override
    public int getHeight() {
        if (!isVisible()) {
            return 0;
        }
        int rowHeight = getDisplayMode() == DisplayMode.SQUISHED ? squishedRowHeight : expandedRowHeight;
        int minHeight = rowHeight * Math.max(1, getNumberOfFeatureLevels());
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

    public int getFeatureWindowSize() {
        return source.getFeatureWindowSize();
    }

    public void setRendererClass(Class rc) {
        try {
            renderer = (FeatureRenderer) rc.newInstance();
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
        if (getDisplayMode() != DisplayMode.COLLAPSED && packedFeaturesMap.size() > 0) {
            int n = 0;
            for (PackedFeatures pf : packedFeaturesMap.values()) {
                //dhmay adding null check.  To my mind this shouldn't be necessary, but we're encountering
                //it intermittently.  Food for future thought
                if (pf != null)
                    n = Math.max(n, pf.getRowCount());
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
            if (scoreType == RegionScoreType.MUTATION_COUNT && this.getTrackType() == TrackType.MUTATION) {
                Iterator<Feature> features = source.getFeatures(chr, start, end);
                int count = 0;
                if (features != null) {
                    while (features.hasNext()) {
                        Feature f = features.next();
                        if (f.getStart() > end) {
                            break;
                        }
                        if (f.getEnd() >= start) {
                            count++;
                        }
                    }
                }
                return count;
            }
        } catch (IOException e) {
            log.error("Error counting features.", e);
        }
        return -Float.MAX_VALUE;
    }


    public FeatureRenderer getRenderer() {
        if (renderer == null) {
            setRendererClass(IGVFeatureRenderer.class);
        }
        return renderer;
    }

    /**
     * Return a string for popup text.
     *
     * @param chr
     * @param position in genomic coordinates
     * @param y        - pixel position in panel coordinates (i.e. not track coordinates)
     * @return
     */
    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {


        if (showFeatures) {

            List<Feature> allFeatures = getAllFeatureAt(position, y, frame);
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
                        buf.append("<br/>--------------<br/>");
                    }

                    IGVFeature igvFeature = (IGVFeature) feature;
                    String vs = igvFeature.getValueString(position, null);
                    buf.append(vs);

                    if (IGV.getInstance().isShowDetailsOnClick()) {
                        // URL
                        String url = getFeatureURL(igvFeature);

                        if (url != null) {
                            buf.append("<br/><a href=\"" + url + "\">" + url + "</a>");
                        }
                    }

                    firstFeature = false;

                    if (n > maxNumber) {
                        buf.append("...");
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


    private String getFeatureURL(IGVFeature igvFeature) {
        String url = igvFeature.getURL();
        if (url == null) {
            String trackURL = getUrl();
            if (trackURL != null && igvFeature.getIdentifier() != null) {
                String encodedID = StringUtils.encodeURL(igvFeature.getIdentifier());
                url = trackURL.replaceAll("\\$\\$", encodedID);
            }
        }
        return url;
    }


    /**
     * Get all features this track contains.
     *
     * @return
     */
    public List<Feature> getFeatures(String chr, int start, int end) {
        List<Feature> features = new ArrayList<Feature>();
        try {
            Iterator<Feature> iter = source.getFeatures(chr, start, end);
            while (iter.hasNext()) {
                features.add(iter.next());
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
    private List<Feature> getAllFeatureAt(double position, int y, ReferenceFrame frame) {

        PackedFeatures<IGVFeature> packedFeatures = packedFeaturesMap.get(frame.getName());

        if (packedFeatures == null) {
            return null;
        }

        List<Feature> feature = null;

        // Determine the level number (for expanded tracks.
        int levelNumber = 0;
        if (levelRects != null) {
            for (int i = 0; i < levelRects.size(); i++) {
                Rectangle r = levelRects.get(i);
                if ((y >= r.y) && (y <= r.getMaxY())) {
                    levelNumber = i;
                    break;
                }
            }
        }

        int nLevels = this.getNumberOfFeatureLevels();
        List<IGVFeature> features = null;
        if ((nLevels > 1) && (levelNumber < nLevels)) {
            features = packedFeatures.getRows().get(levelNumber).getFeatures();
        } else {
            features = packedFeatures.getFeatures();
        }
        if (features != null) {

            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            double minWidth = MINIMUM_FEATURE_SPACING * bpPerPixel;
            // The maximum length of all features in this collection. Used to insure we consider all features that
            // might overlap the position (feature are sorted by start position, but length is variable)
            int maxFeatureLength = packedFeatures.getMaxFeatureLength();
            feature = FeatureUtils.getAllFeaturesAt(position, maxFeatureLength, minWidth, features);
        }
        return feature;
    }

    /**
     * Determine which row the user clicked in and return the appropriate feature
     *
     * @param chr
     * @param position
     * @param y
     * @param frame
     * @return
     */
    public Feature getFeatureAt(String chr, double position, int y, ReferenceFrame frame) {
        // Determine the level number (for expanded tracks.
        int featureRow = 0;
        if (levelRects != null) {
            for (int i = 0; i < levelRects.size(); i++) {
                Rectangle r = levelRects.get(i);
                if ((y >= r.y) && (y <= r.getMaxY())) {
                    featureRow = i;
                    break;
                }
            }
        }
        return getFeatureAtPositionInFeatureRow(position, featureRow, frame);
    }

    /**
     * Knowing the feature row, figure out which feature is at position position. If not expanded,
     * featureRow is ignored
     *
     * @param position
     * @param featureRow
     * @param frame
     * @return
     */
    public Feature getFeatureAtPositionInFeatureRow(double position, int featureRow, ReferenceFrame frame) {

        PackedFeatures<IGVFeature> packedFeatures = packedFeaturesMap.get(frame.getName());

        if (packedFeatures == null) {
            return null;
        }

        Feature feature = null;

        int nLevels = this.getNumberOfFeatureLevels();
        List<IGVFeature> features = null;
        if ((nLevels > 1) && (featureRow < nLevels)) {
            features = packedFeatures.getRows().get(featureRow).getFeatures();
        } else {
            features = packedFeatures.getFeatures();
        }
        if (features != null) {

            // give a +/- 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            int minWidth = (int) (2 * bpPerPixel);
            feature = FeatureUtils.getFeatureAt(position, minWidth, features);
        }
        return feature;
    }


    public WindowFunction getWindowFunction() {
        return WindowFunction.count;
    }

    @Override
    public boolean handleDataClick(TrackClickEvent te) {

        MouseEvent e = te.getMouseEvent();

        //dhmay adding selection of an expanded feature row
        if (getDisplayMode() != DisplayMode.COLLAPSED) {
            if (levelRects != null) {
                for (int i = 0; i < levelRects.size(); i++) {
                    Rectangle rect = levelRects.get(i);
                    if (rect.contains(e.getPoint())) {
                        if (i == selectedFeatureRowIndex)
                            setSelectedFeatureRowIndex(FeatureTrack.NO_FEATURE_ROW_SELECTED);
                        else {
                            //make this track selected
                            setSelected(true);
                            //select the appropriate row
                            setSelectedFeatureRowIndex(i);
                        }
                        IGV.getInstance().doRefresh();
                        break;
                    }
                }
            }
        }


        //dhmay adding for feature selection
        selectedFeature = null;

        Feature f = getFeatureAtMousePosition(te);
        if (f != null && f instanceof IGVFeature) {
            IGVFeature igvFeature = (IGVFeature) f;
            //if nothing already selected, select this feature
            if (selectedFeature == null)
                selectedFeature = igvFeature;
                //If something already selected, then if it's the same as this feature, deselect, otherwise, select
                //this feature.
                //todo: contains() might not do everything I want it to.
            else if (igvFeature.contains(selectedFeature) && (selectedFeature.contains(igvFeature)))
                selectedFeature = null;
            else selectedFeature = igvFeature;

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
            double location = referenceFrame.getChromosomePosition(e.getX());
            Feature f = getFeatureAt(referenceFrame.getChrName(), location, e.getY(), referenceFrame);
            return f;
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
        super.setDisplayMode(mode);
    }

    @Override
    public void preload(RenderContext context) {
        ReferenceFrame frame = context.getReferenceFrame();
        PackedFeatures packedFeatures = packedFeaturesMap.get(frame.getName());
        String chr = context.getChr();
        int start = (int) context.getOrigin();
        int end = (int) context.getEndLocation();
        if (packedFeatures == null || !packedFeatures.containsInterval(chr, start, end)) {
            loadFeatures(frame.getChrName(), (int) frame.getOrigin(), (int) frame.getEnd(), context);
        }
    }


    public void render(RenderContext context, Rectangle rect) {
        Rectangle renderRect = new Rectangle(rect);
        renderRect.y = renderRect.y + margin;
        renderRect.height -= margin;


        showFeatures = isShowFeatures(context);
        if (showFeatures) {
            if (lastFeatureMode != null) {
                super.setDisplayMode(lastFeatureMode);
                lastFeatureMode = null;
            }
            renderFeatures(context, renderRect);
        } else {
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
            borderGraphics.drawLine(rect.x, rect.y + rect.height, rect.x + rect.width, rect.y + rect.height);
        }

    }

    protected boolean isShowFeatures(RenderContext context) {
        double windowSize = context.getEndLocation() - context.getOrigin();
        int vw = getVisibilityWindow();
        return (vw <= 0 && !context.getChr().equals(Globals.CHR_ALL) ||
                windowSize <= vw && !context.getChr().equals(Globals.CHR_ALL));
    }

    protected void renderCoverage(RenderContext context, Rectangle inputRect) {
        if (source == null) {
            return;
        }

        List<LocusScore> scores = source.getCoverageScores(context.getChr(), (int) context.getOrigin(),
                (int) context.getEndLocation(), context.getZoom());
        if (scores == null) {
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            Rectangle textRect = new Rectangle(inputRect);

            // Keep text near the top of the track rectangle
            textRect.height = Math.min(inputRect.height, 20);
            String message = context.getChr().equals(Globals.CHR_ALL) ? "Zoom in to see features." :
                    "Zoom in to see features, or right-click to increase Feature Visibility Window.";
            GraphicUtils.drawCenteredText(message, textRect, g);

        } else {
            float max = getMaxEstimate(scores);
            ContinuousColorScale cs = getColorScale();
            if (cs != null) {
                cs.setPosEnd(max);
            }
            setDataRange(new DataRange(0, 0, max));
            coverageRenderer.render(scores, context, inputRect, this);
        }
    }

    private float getMaxEstimate(List<LocusScore> scores) {
        float max = 0;
        int n = Math.min(200, scores.size());
        for (int i = 0; i < n; i++) {
            max = Math.max(max, scores.get(i).getScore());
        }
        return max;
    }

    // Render features in the given input rectangle.

    protected void renderFeatures(RenderContext context, Rectangle inputRect) {

        if (featuresLoading || fatalLoadError) {
            return;
        }

        if (log.isDebugEnabled()) {
            log.debug("renderFeatures: " + getName());
        }

        String chr = context.getChr();
        int start = (int) context.getOrigin();
        int end = (int) context.getEndLocation() + 1;

        PackedFeatures packedFeatures = packedFeaturesMap.get(context.getReferenceFrame().getName());

        if (packedFeatures == null || !packedFeatures.containsInterval(chr, start, end)) {
            loadFeatures(chr, start, end, context);
            if (!IGV.hasInstance() || !IGV.getInstance().isExportingSnapshot()) {
                // DONT CALL REPAINT HERE!!! FEATURES ARE LOADING ASYNCHRONOUSLY, REPAINT CALLED WHEN LOADING IS DONE
                return;
            }
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
                    IGV.getInstance().removeTracks(tmp);
                    IGV.getInstance().doRefresh();
                } else {
                    fatalLoadError = false;
                }
            }
        }


    }

    protected void renderFeatureImpl(RenderContext context, Rectangle inputRect, PackedFeatures packedFeatures) {


        FeatureRenderer renderer = getRenderer();
        if (getDisplayMode() != DisplayMode.COLLAPSED) {
            List<PackedFeatures.FeatureRow> rows = packedFeatures.getRows();
            if (rows != null && rows.size() > 0) {

                int nLevels = rows.size();
                synchronized (levelRects) {

                    levelRects.clear();

                    // Divide rectangle into equal height levels
                    double h = inputRect.getHeight() / nLevels;
                    Rectangle rect = new Rectangle(inputRect.x, inputRect.y, inputRect.width, (int) h);
                    int i = 0;

                    renderer.reset();
                    for (PackedFeatures.FeatureRow row : rows) {
                        levelRects.add(new Rectangle(rect));
                        renderer.render(row.features, context, levelRects.get(i), this);
                        if (selectedFeatureRowIndex == i) {
                            Graphics2D fontGraphics = context.getGraphic2DForColor(SELECTED_FEATURE_ROW_COLOR);
                            fontGraphics.fillRect(rect.x, rect.y, rect.width, rect.height);
                        }
                        rect.y += h;
                        i++;
                    }
                }
            }
        } else {
            List<IGVFeature> features = packedFeatures.getFeatures();
            if (features != null) {
                renderer.render(features, context, inputRect, this);
            }
        }
    }


    /**
     * Loads and segregates features into rows such that they do not overlap.  Loading is done in a background
     * thread.
     *
     * @param chr
     * @param start
     * @param end
     */
    protected synchronized void loadFeatures(final String chr, final int start, final int end, final RenderContext context) {

        featuresLoading = true;

        // TODO -- improve or remove the need for this test.  We know that FeatureCollectionSource has all the data
        // in memory, and can by run synchronously
        boolean aSync = !(source instanceof FeatureCollectionSource);


        NamedRunnable runnable = new NamedRunnable() {
            public void run() {
                try {
                    featuresLoading = true;

                    int maxEnd = end;
                    Genome genome = GenomeManager.getInstance().getCurrentGenome();
                    if (genome != null) {
                        Chromosome c = genome.getChromosome(chr);
                        if (c != null) maxEnd = Math.max(c.getLength(), end);
                    }
                    int delta = (end - start) / 2;
                    int expandedStart = start - delta;
                    int expandedEnd = end + delta;


                    // TODO -- implement source to return iterators
                    Iterator<Feature> iter = source.getFeatures(chr, expandedStart, expandedEnd);
                    if (iter == null) {
                        PackedFeatures pf = new PackedFeatures(chr, expandedStart, expandedEnd);
                        packedFeaturesMap.put(context.getReferenceFrame().getName(), pf);
                    } else {
                        //dhmay putting a switch in for different packing behavior in splice junction tracks.
                        //This should probably be switched somewhere else, but that would require a big refactor.
                        PackedFeatures pf = new PackedFeatures(chr, expandedStart, expandedEnd, iter, getName());
                        packedFeaturesMap.put(context.getReferenceFrame().getName(), pf);
                    }


                    if (IGV.hasInstance()) {
                        // TODO -- WHY IS THIS HERE????
                        IGV.getInstance().layoutMainPanel();
                    }
                    if (context.getPanel() != null) context.getPanel().repaint();
                } catch (Exception e) {
                    // Mark the interval with an empty feature list to prevent an endless loop of load
                    // attempts.
                    PackedFeatures pf = new PackedFeatures(chr, start, end);
                    packedFeaturesMap.put(context.getReferenceFrame().getName(), pf);
                    String msg = "Error loading features for interval: " + chr + ":" + start + "-" + end + " <br>" + e.toString();
                    MessageUtils.showMessage(msg);
                    log.error(msg, e);
                } finally {
                    featuresLoading = false;
                }
            }

            public String getName() {
                return "Load features: " + FeatureTrack.this.getName();
            }
        };

        if (aSync) {
            LongRunningTask.submit(runnable);
        } else {
            runnable.run();
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

        if (packedFeatures != null && packedFeatures.containsInterval(chr, (int) center - 1, (int) center + 1)) {
            if (packedFeatures.getFeatures().size() > 0 && canScroll) {
                f = (forward ? FeatureUtils.getFeatureAfter(center + 1, packedFeatures.getFeatures()) :
                        FeatureUtils.getFeatureBefore(center - 1, packedFeatures.getFeatures()));
            }

            if (f == null) {
                int binSize = source.getFeatureWindowSize();

                final Genome genome = GenomeManager.getInstance().getCurrentGenome();
                if (forward == true) {
                    // Forward
                    int nextStart = packedFeatures.getEnd();
                    String nextChr = chr;
                    while (nextChr != null) {
                        int chrLength = genome.getChromosome(nextChr).getLength();
                        while (nextStart < chrLength) {
                            int nextEnd = binSize > 0 ? nextStart + source.getFeatureWindowSize() : chrLength;
                            Iterator<Feature> iter = source.getFeatures(nextChr, nextStart, nextEnd);
                            // The check on position should not be neccessary, but not all implementations of getFeatures
                            // obey the contract to return features only in the interval.
                            if (iter != null) {
                                while (iter.hasNext()) {
                                    Feature feat = iter.next();
                                    if (feat.getStart() > nextStart) {
                                        return feat;
                                    }
                                }
                            }
                            nextStart = nextEnd;
                        }
                        nextChr = genome.getNextChrName(nextChr);
                        nextStart = 0;
                    }
                } else {
                    // Reverse
                    int nextEnd = packedFeatures.getStart();
                    String nextChr = chr;
                    while (nextChr != null) {
                        while (nextEnd > 0) {
                            int nextStart = binSize > 0 ? Math.max(0, nextEnd - source.getFeatureWindowSize()) : 0;
                            Iterator<Feature> iter = source.getFeatures(nextChr, nextStart, nextEnd);
                            if (iter != null && iter.hasNext()) {
                                // The check on position should not be neccessary, but not all implementations of getFeatures
                                // obey the contract to return features only in the interval.

                                Feature prevFeature = null;
                                while (iter.hasNext()) {
                                    Feature feat = iter.next();
                                    if (feat.getStart() < nextEnd) {
                                        prevFeature = feat;
                                    }
                                }
                                if (prevFeature != null) {
                                    return prevFeature;
                                }
                            }
                            nextEnd = nextStart;
                        }
                        nextChr = genome.getPrevChrName(nextChr);
                        if (nextChr != null) {
                            nextEnd = genome.getChromosome(nextChr).getLength();
                        }
                    }
                }
            }
        }

        return f;
    }

    public void setVisibilityWindow(int windowSize) {
        super.setVisibilityWindow(windowSize);
        packedFeaturesMap.clear();
        source.setFeatureWindowSize(visibilityWindow);
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

    @Override
    public void restorePersistentState(Map<String, String> attributes) {
        super.restorePersistentState(attributes);    //To change body of overridden methods use File | Settings | File Templates.


    }

    @Override
    public Map<String, String> getPersistentState() {
        Map<String, String> stateMap = super.getPersistentState();
        return stateMap;

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
}

