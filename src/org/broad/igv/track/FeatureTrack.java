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
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.renderer.*;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.tribble.Feature;

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

    final String FEATURE_VISIBILITY_WINDOW = "featureVisibilityWindow";

    static int maxLevels = 200;

    private boolean expanded;
    private List<Rectangle> featureRects = new ArrayList();

    // TODO -- this is a memory leak
    protected LRUCache<String, PackedFeatures> packedFeaturesMap = new LRUCache(200);

    private FeatureRenderer renderer;
    private DataRenderer coverageRenderer;

    // true == features,  false =  coverage
    private boolean showFeatures = true;

    protected FeatureSource source;
    public static final int MINIMUM_FEATURE_SPACING = 1;


    private boolean featuresLoading = false;

    private Rectangle EXPAND_BUTTON_RECT = new Rectangle();
    private static final int EXPAND_ICON_BUFFER_WIDTH = 17;
    private static final int EXPAND_ICON_BUFFER_HEIGHT = 17;
    public static final int MARGIN = 5;

    //track which row of the expanded track is selected by the user.
    //Selection goes away if tracks are collpased
    public static final int NO_FEATURE_ROW_SELECTED = -1;
    private static final Color SELECTED_FEATURE_ROW_COLOR = new Color(50, 170, 50, 30);
    int selectedFeatureRowIndex = NO_FEATURE_ROW_SELECTED;


    public FeatureTrack(String id, FeatureSource source) {
        super(id);
        init(source);

    }


    public FeatureTrack(ResourceLocator locator, FeatureSource source) {
        super(locator);
        init(source);
        this.getMinimumHeight();
    }

    public FeatureTrack(ResourceLocator locator) {
        super(locator);
    }

    private void init(FeatureSource source) {
        this.expanded = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.EXPAND_FEAUTRE_TRACKS);
        this.source = source;
        setMinimumHeight(10);
        setColor(Color.blue.darker());
        renderer = new IGVFeatureRenderer();
        coverageRenderer = new HeatmapRenderer();
        if (source.getFeatureWindowSize() > 0) {
            visibilityWindow = source.getFeatureWindowSize();
        }
    }

    public void setRendererClass(Class rc) {
        try {
            renderer = (FeatureRenderer) rc.newInstance();
        } catch (Exception ex) {
            log.error("Error instatiating renderer ", ex);
        }
    }

    @Override
    public void setTrackProperties(TrackProperties trackProperties) {
        super.setTrackProperties(trackProperties);

        if (trackProperties.getFeatureVisibilityWindow() >= 0) {
            setVisibilityWindow(trackProperties.getFeatureVisibilityWindow());
        }

    }

    @Override
    public int getHeight() {

        if (false == isVisible()) {
            return 0;
        }
        return super.getHeight() * Math.max(1, getNumberOfFeatureLevels());
    }

    public int getNumberOfFeatureLevels() {
        if (expanded && packedFeaturesMap.size() > 0) {
            int n = 0;
            for (PackedFeatures pf : packedFeaturesMap.values()) {
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
     * @param y - pixel position in panel coordinates (i.e. not track coordinates)
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
                // give a 2 pixel window, otherwise very narrow features will be missed.
                double bpPerPixel = frame.getScale();
                double minWidth = MINIMUM_FEATURE_SPACING * bpPerPixel;    /* * */
                LocusScore score = (LocusScore) FeatureUtils.getFeatureAt(position, minWidth, scores);
                return score == null ? "" : "Mean count: " + score.getScore();
            }

        }
    }

    /**
     *
     * @param position in genomic coordinates
     * @param y  pixel location in panel coordinates.  // TODO offset by track origin before getting here?
     * @param frame
     * @return
     */
    protected List<Feature> getAllFeatureAt(String chr, double position, int y, ReferenceFrame frame) {

        PackedFeatures packedFeatures = packedFeaturesMap.get(frame.getName());

        if (packedFeatures == null) {
            return null;
        }

        List<Feature> feature = null;
        // Determine the level number (for expanded tracks.
        int levelNumber = 0;
        if (featureRects != null) {
            for (int i = 0; i < featureRects.size(); i++) {
                Rectangle r = featureRects.get(i);
                if ((y >= r.y) && (y <= r.getMaxY())) {
                    levelNumber = i;
                    break;
                }
            }
        }

        int nLevels = this.getNumberOfFeatureLevels();
        List<Feature> features = null;
        if ((nLevels > 1) && (levelNumber < nLevels)) {
            features = packedFeatures.getRows().get(levelNumber).features;
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
     * @param chr
     * @param position
     * @param y
     * @param frame
     * @return
     */
    public Feature getFeatureAt(String chr, double position, int y, ReferenceFrame frame) {
        // Determine the level number (for expanded tracks.
        int featureRow = 0;
        if (featureRects != null) {
            for (int i = 0; i < featureRects.size(); i++) {
                Rectangle r = featureRects.get(i);
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
     * @param chr
     * @param position
     * @param featureRow
     * @param frame
     * @return
     */
    public Feature getFeatureAtPositionInFeatureRow(String chr, double position, int featureRow, 
                                                    ReferenceFrame frame) {

        PackedFeatures packedFeatures = packedFeaturesMap.get(frame.getName());

        if (packedFeatures == null) {
            return null;
        }

        Feature feature = null;

        int nLevels = this.getNumberOfFeatureLevels();
        List<Feature> features = null;
        if ((nLevels > 1) && (featureRow < nLevels)) {
            features = packedFeatures.getRows().get(featureRow).features;
        } else {
            features = packedFeatures.getFeatures();
        }
        if (features != null) {

            // give a 2 pixel window, otherwise very narrow features will be missed.
            double bpPerPixel = frame.getScale();
            double minWidth = MINIMUM_FEATURE_SPACING * bpPerPixel;
            feature = FeatureUtils.getFeatureAt(position, minWidth, features, true);
        }
        return feature;
    }

    protected Feature getFeatureClosest(double position, int y, ReferenceFrame frame) {

        PackedFeatures packedFeatures = packedFeaturesMap.get(frame.getName());

        if (packedFeatures == null) {
            return null;
        }

        Feature feature = null;
        // Determine the level number (for expanded tracks.
        int levelNumber = 0;
        if (featureRects != null) {
            for (int i = 0; i < featureRects.size(); i++) {
                Rectangle r = featureRects.get(i);
                if ((y >= r.y) && (y <= r.getMaxY())) {
                    levelNumber = i;
                    break;
                }
            }
        }

        int nLevels = this.getNumberOfFeatureLevels();
        List<Feature> features = null;
        if ((nLevels > 1) && (levelNumber < nLevels)) {
            features = packedFeatures.getRows().get(levelNumber).features;
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


        if (e.getX() < EXPAND_ICON_BUFFER_WIDTH) {
            if (EXPAND_BUTTON_RECT.contains(e.getPoint())) {
                setExpanded(!expanded);
                IGVMainFrame.getInstance().doRefresh();
            }
            return true;
        }

        //dhmay adding selection of an expanded feature row
        if (expanded)
        {
            if (featureRects != null) {
                for (int i=0; i<featureRects.size(); i++)
                {
                    Rectangle rect = featureRects.get(i);
                    if (rect.contains(e.getPoint()))
                    {
                        if (i == selectedFeatureRowIndex)
                            setSelectedFeatureRowIndex(FeatureTrack.NO_FEATURE_ROW_SELECTED);
                        else
                            setSelectedFeatureRowIndex(i);
                        IGVMainFrame.getInstance().doRefresh();
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


    @Override
    public boolean isExpanded() {
        return expanded;
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

    public void render(RenderContext context, Rectangle rect) {
        Rectangle renderRect = new Rectangle(rect);
        renderRect.y = renderRect.y + MARGIN;
        renderRect.height -= MARGIN;


        double windowSize = context.getEndLocation() - context.getOrigin();

        int vw = getVisibilityWindow();
        showFeatures = (vw <= 0 && !context.getChr().equals(Globals.CHR_ALL) ||
                windowSize <= vw && !context.getChr().equals(Globals.CHR_ALL));
        if (showFeatures) {
            renderFeatures(context, renderRect);
        } else {
            renderCoverage(context, renderRect);
        }
        boolean showIcon = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_EXPAND_ICON);
        if (showIcon && !IGVMainFrame.getInstance().isExportingSnapshot() && !FrameManager.isGeneListMode()) {
            renderExpandTool(context, rect);
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
            coverageRenderer.render(this, scores, context, inputRect);
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


    private int[] p1 = new int[3];
    private int[] p2 = new int[3];

    private void renderExpandTool(RenderContext contect, Rectangle rect) {

        PackedFeatures packedFeatures = packedFeaturesMap.get(contect.getReferenceFrame().getName());

        if (packedFeatures == null || packedFeatures.getRowCount() <= 1) {
            return;
        }

        Graphics2D g2d = contect.getGraphic2DForColor(Color.DARK_GRAY);
        int levelHeight = getHeight() / (this.getNumberOfFeatureLevels() + 1);

        g2d.clearRect(rect.x, rect.y, EXPAND_ICON_BUFFER_WIDTH, levelHeight);

        EXPAND_BUTTON_RECT.x = rect.x + 3;
        EXPAND_BUTTON_RECT.y = rect.y + MARGIN + 4;
        EXPAND_BUTTON_RECT.width = 10;
        EXPAND_BUTTON_RECT.height = 10;

        if (expanded) {
            p1[0] = EXPAND_BUTTON_RECT.x;
            p1[1] = EXPAND_BUTTON_RECT.x + 8;
            p1[2] = EXPAND_BUTTON_RECT.x + 4;
            p2[0] = EXPAND_BUTTON_RECT.y;
            p2[1] = EXPAND_BUTTON_RECT.y;
            p2[2] = EXPAND_BUTTON_RECT.y + 8;
            g2d.fillPolygon(p1, p2, 3);

        } else {
            p1[0] = EXPAND_BUTTON_RECT.x;
            p1[1] = EXPAND_BUTTON_RECT.x + 8;
            p1[2] = EXPAND_BUTTON_RECT.x;
            p2[0] = EXPAND_BUTTON_RECT.y;
            p2[1] = EXPAND_BUTTON_RECT.y + 4;
            p2[2] = EXPAND_BUTTON_RECT.y + 8;
            g2d.fillPolygon(p1, p2, 3);
        }

    }


    // Render features in the given input rectangle.

    private void renderFeatures(RenderContext context, Rectangle inputRect) {

        if (featuresLoading) {
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

            if (!IGVMainFrame.getInstance().isExportingSnapshot()) {
                return;
            }


        }


        renderFeatureImpl(context, inputRect, packedFeatures);


    }

    protected void renderFeatureImpl(RenderContext context, Rectangle inputRect, PackedFeatures packedFeatures) {
        if (expanded) {
            List<FeatureRow> rows = packedFeatures.getRows();
            if (rows != null && rows.size() > 0) {

                int nLevels = rows.size();
                synchronized (featureRects) {

                    featureRects.clear();

                    // Divide rectangle into equal height levels
                    double h = inputRect.getHeight() / nLevels;
                    Rectangle rect = new Rectangle(inputRect.x, inputRect.y, inputRect.width, (int) h);
                    int i=0;
                    for (FeatureRow row : rows) {
                        featureRects.add(new Rectangle(rect));
                        getRenderer().renderFeatures(row.features, context, rect, this);
                        if (selectedFeatureRowIndex == i)
                        {
                            Graphics2D fontGraphics =
                                    (Graphics2D) context.getGraphic2DForColor(SELECTED_FEATURE_ROW_COLOR).create();
                            fontGraphics.fillRect(rect.x, rect.y, rect.width, rect.height);
                        }
                        rect.y += h;
                        i++;
                    }
                }
            }
        } else {
            List<Feature> features = packedFeatures.getFeatures();
            if (features != null) {
                getRenderer().renderFeatures(features, context, inputRect, this);
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
    private synchronized void loadFeatures(final String chr, final int start, final int end, final RenderContext context) {

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
                    int expandedStart = Math.max(0, start - delta);
                    int expandedEnd = Math.min(maxEnd, end + delta);


                    // TODO -- implement source to return iterators
                    Iterator<Feature> iter = source.getFeatures(chr, start, end);
                    if (iter == null) {
                        packedFeaturesMap.put(context.getReferenceFrame().getName(), new PackedFeatures(chr, expandedStart, expandedEnd));
                    } else {
                        packedFeaturesMap.put(context.getReferenceFrame().getName(), new PackedFeatures(chr, expandedStart, expandedEnd, iter, getName()));
                    }

                    IGVMainFrame.getInstance().layoutMainPanel();
                    if (context.getPanel() != null) context.getPanel().repaint();
                } catch (Throwable e) {
                    // Mark the interval with an empty feature list to prevent an endless loop of load
                    // attempts.
                    packedFeaturesMap.put(context.getReferenceFrame().getName(), new PackedFeatures(chr, start, end));
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


    @Override
    public void setExpanded(boolean value) {
        expanded = value;
        //dhmay adding: if collapsing rows, remove row selection if any
        if (!expanded)
            setSelectedFeatureRowIndex(NO_FEATURE_ROW_SELECTED);
        IGVMainFrame.getInstance().layoutMainPanel();   // TODO -- this is for the scrollbar.  Can some event be fired?
    }

    @Override
    public void setHeight(int newHeight) {

        int levelCount = this.getNumberOfFeatureLevels();
        super.setHeight(Math.max(getMinimumHeight(), newHeight / levelCount));
    }


    public void setStatType(WindowFunction type) {
    }

    /**
     * Method description
     *
     * @param zoom
     */
    public void setZoom(int zoom) {
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

    @Override
    public void restorePersistentState(Map<String, String> attributes) {
        super.restorePersistentState(attributes);    //To change body of overridden methods use File | Settings | File Templates.

        String fvw = attributes.get(FEATURE_VISIBILITY_WINDOW);
        if (fvw != null) {
            try {
                visibilityWindow = Integer.parseInt(fvw);
            } catch (NumberFormatException e) {
                log.error("Error restoring featureVisibilityWindow: " + fvw);
            }
        }

    }

    @Override
    public Map<String, String> getPersistentState() {
        Map<String, String> stateMap = super.getPersistentState();
        stateMap.put(FEATURE_VISIBILITY_WINDOW, String.valueOf(visibilityWindow));
        return stateMap;

    }

    public int getSelectedFeatureRowIndex() {
        return selectedFeatureRowIndex;
    }

    public void setSelectedFeatureRowIndex(int selectedFeatureRowIndex) {
        this.selectedFeatureRowIndex = selectedFeatureRowIndex;
    }

    //public IGVFeature nextFeature(String chr, double position, boolean forward) {
//
//    return source.nextFeature(chr, position, forward);
//}


    static class FeatureRow {
        int start;
        int end;
        List<Feature> features;

        FeatureRow() {
            this.features = new ArrayList(100);
        }

        void addFeature(Feature feature) {
            if (features.isEmpty()) {
                this.start = feature.getStart();
            }
            features.add(feature);
            end = feature.getEnd();
        }
    }


}

