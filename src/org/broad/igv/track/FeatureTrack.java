/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.renderer.*;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.tribble.Feature;
import org.broad.tribble.TribbleException;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.net.URLEncoder;
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
    int selectedFeatureRowIndex = NO_FEATURE_ROW_SELECTED;

    int margin = DEFAULT_MARGIN;

    private static boolean drawBorder = true;

    /**
     * Does not initialize with the featuresource
     *
     * @param id
     * @param name
     */
    public FeatureTrack(String id, String name) {
        super(id, name);
        setSortable(false);
    }

    public FeatureTrack(String id, String name, FeatureSource source) {
        super(id, name);
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
        this.getMinimumHeight();
        setSortable(false);
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
    public void chromosomeChanged(String chrName) {
        if (chrName.equals(Globals.CHR_ALL)) {
            showFeatures = false;
        } else if (showFeatures == false) {
            // Currently showing coverage data.
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


    @Override
    public void setHeight(int newHeight) {
        super.setHeight(newHeight);
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

    }


    public void setWindowFunction(WindowFunction type) {
        // Ignored for feature tracks
    }

    public int getNumberOfFeatureLevels() {
        if (getDisplayMode() != DisplayMode.COLLAPSED && packedFeaturesMap.size() > 0) {
            int n = 0;
            for (PackedFeatures pf : packedFeaturesMap.values()) {
                //dhmay adding null check.  To my mind this shouldn't be necessary, but we're encountering
                //it intermittently.  Food for future thought
                if (pf != null)
                    n += pf.getRowCount();
            }
            return n;
        }
        return 1;
    }

    /**
     * Return a score over the interval.  This is required by the track interface to support sorting.
     */
    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType scoreType, ReferenceFrame frame) {

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

            List<Feature> allFeatures = getAllFeatureAt(chr, position, y, frame);
            if (allFeatures == null) {
                return "";
            }

            StringBuffer buf = new StringBuffer();
            boolean firstFeature = true;
            int maxNumber = 10;
            int n = 1;
            for (Feature feature : allFeatures) {
                if (feature != null && feature instanceof IGVFeature) {
                    IGVFeature igvFeature = (IGVFeature) feature;
                    String vs = igvFeature.getValueString(position, null);
                    if (!firstFeature) {
                        buf.append("<br>--------------<br>");
                    }
                    buf.append(vs);
                    firstFeature = false;

                    if (n > maxNumber) {
                        buf.append("...");
                        break;
                    }
                }
                n++;
            }
            if (!firstFeature) buf.append("<br>--------------<br>");

            return buf.toString();
        } else {
            int zoom = Math.max(0, frame.getZoom());
            List<LocusScore> scores = source.getCoverageScores(chr, (int) position - 10, (int) position + 10, zoom);

            if (scores == null) {
                return "";
            } else {
                // give a +/- 2 pixel buffer, otherwise very narrow features will be missed.
                double bpPerPixel = frame.getScale();
                int minWidth = (int) (2 * bpPerPixel);    /* * */
                LocusScore score = (LocusScore) FeatureUtils.getFeatureAt(position, minWidth, scores);
                return score == null ? "" : "Mean count: " + score.getScore();
            }

        }
    }

    /**
     * @param position in genomic coordinates
     * @param y        pixel location in panel coordinates.  // TODO offset by track origin before getting here?
     * @param frame
     * @return
     */
    protected List<Feature> getAllFeatureAt(String chr, double position, int y, ReferenceFrame frame) {

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
            feature = FeatureUtils.getAllFeaturesAt(position, maxFeatureLength, minWidth, features, true);
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
        return getFeatureAtPositionInFeatureRow(chr, position, featureRow, frame);
    }

    /**
     * Knowing the feature row, figure out which feature is at position position. If not expanded,
     * featureRow is ignored
     *
     * @param chr
     * @param position
     * @param featureRow
     * @param frame
     * @return
     */
    public Feature getFeatureAtPositionInFeatureRow(String chr, double position, int featureRow,
                                                    ReferenceFrame frame) {

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
            feature = FeatureUtils.getFeatureAt(position, minWidth, features, true);
        }
        return feature;
    }

    protected Feature getFeatureClosest(double position, int y, ReferenceFrame frame) {

        PackedFeatures<IGVFeature> packedFeatures = packedFeaturesMap.get(frame.getName());

        if (packedFeatures == null) {
            return null;
        }

        Feature feature = null;
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
            feature = FeatureUtils.getFeatureClosest(position, features);
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


        Feature f = getFeatureAtMousePosition(te);
        if (f != null && f instanceof IGVFeature) {
            IGVFeature igvFeature = (IGVFeature) f;
            String url = igvFeature.getURL();
            if (url == null) {
                String trackURL = getUrl();
                if (trackURL != null && igvFeature.getIdentifier() != null) {
                    String encodedID = URLEncoder.encode(igvFeature.getIdentifier());
                    url = trackURL.replaceAll("\\$\\$", encodedID);
                }
            }

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

        return false;
    }

    public Feature getFeatureAtMousePosition(TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();
        final ReferenceFrame referenceFrame = te.getFrame();
        if (referenceFrame != null) {
            double location = referenceFrame.getChromosomePosition(e.getX());
            double displayLocation = location + 1;
            Feature f = getFeatureAt(referenceFrame.getChrName(), displayLocation, e.getY(), referenceFrame);
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
        super.setDisplayMode(mode);    //To change body of overridden methods use File | Settings | File Templates.
    }

    public void render(RenderContext context, Rectangle rect) {
        Rectangle renderRect = new Rectangle(rect);
        renderRect.y = renderRect.y + margin;
        renderRect.height -= margin;


        double windowSize = context.getEndLocation() - context.getOrigin();

        int vw = getVisibilityWindow();
        showFeatures = (vw <= 0 && !context.getChr().equals(Globals.CHR_ALL) ||
                windowSize <= vw && !context.getChr().equals(Globals.CHR_ALL));
        if (showFeatures) {
            if (lastFeatureMode != null) {
                super.setDisplayMode(lastFeatureMode);
                lastFeatureMode = null;
            }
            renderFeatures(context, renderRect);
        } else {
            if (getDisplayMode() != DisplayMode.COLLAPSED) {
                lastFeatureMode = getDisplayMode();
                super.setDisplayMode(DisplayMode.COLLAPSED);
            }
            renderCoverage(context, renderRect);
        }

        if (FeatureTrack.drawBorder) {
            Graphics2D borderGraphics = context.getGraphic2DForColor(UIConstants.TRACK_BORDER_GRAY);
            borderGraphics.draw(rect);
        }

    }

    private void renderCoverage(RenderContext context, Rectangle inputRect) {
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
            featuresLoading = true;
            loadFeatures(chr, start, end, context);
            if (!IGV.getInstance().isExportingSnapshot()) {
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
                    IGV.getInstance().getTrackManager().removeTracks(tmp);
                    IGV.getInstance().doRefresh();
                } else {
                    fatalLoadError = false;
                }
            }
        }


    }

    protected void renderFeatureImpl(RenderContext context, Rectangle inputRect, PackedFeatures packedFeatures) {
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
                    for (PackedFeatures.FeatureRow row : rows) {
                        levelRects.add(new Rectangle(rect));
                        getRenderer().render(row.features, context, rect, this);
                        if (selectedFeatureRowIndex == i) {
                            Graphics2D fontGraphics =
                                    (Graphics2D) context.getGraphic2DForColor(SELECTED_FEATURE_ROW_COLOR);
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
                getRenderer().render(features, context, inputRect, this);
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

        // TODO -- improve or remove the need for this test.  We know that FeatureCollectionSource has all the data
        // in memory, and can by run synchronously
        boolean aSync = !(source instanceof FeatureCollectionSource);


        NamedRunnable runnable = new NamedRunnable() {
            public void run() {
                try {
                    featuresLoading = true;

                    int maxEnd = end;
                    Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
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
                        PackedFeatures pf = null;
                        if (getRenderer() instanceof SpliceJunctionRenderer)
                            pf = new PackedFeaturesSpliceJunctions(chr, expandedStart, expandedEnd, iter, getName());
                        else
                            pf = new PackedFeatures(chr, expandedStart, expandedEnd, iter, getName());
                        packedFeaturesMap.put(context.getReferenceFrame().getName(), pf);
                    }

                    IGV.getInstance().layoutMainPanel();
                    if (context.getPanel() != null) context.getPanel().repaint();
                } catch (Throwable e) {
                    // Mark the interval with an empty feature list to prevent an endless loop of load
                    // attempts.
                    PackedFeatures pf = new PackedFeatures(chr, start, end);
                    packedFeaturesMap.put(context.getReferenceFrame().getName(), pf);
                    String msg = "Error loading features for interval: " + chr + ":" + start + "-" + end + " <br>" + e.toString();
                    MessageUtils.showMessage(msg);
                    log.error(msg, e);
                }
                finally {
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

                final Genome genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
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
}

