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

package org.broad.igv.sam;


import org.broad.igv.Globals;
import org.broad.igv.event.AlignmentTrackEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.event.IGVEventObserver;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.ChromosomeNameComparator;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.jbrowse.CircularViewUtilities;
import org.broad.igv.lists.GeneList;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.sashimi.SashimiPlot;
import org.broad.igv.session.Persistable;
import org.broad.igv.session.Session;
import org.broad.igv.tools.PFMExporter;
import org.broad.igv.track.*;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.InsertSizeSettingsDialog;
import org.broad.igv.ui.color.ColorTable;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.Pair;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.util.collections.CollUtils;
import org.broad.igv.util.extview.ExtendViewClient;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.util.List;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */

public class AlignmentTrack extends AbstractTrack implements IGVEventObserver {

    // Alignment colors
    static Color DEFAULT_ALIGNMENT_COLOR = new Color(185, 185, 185); //200, 200, 200);

    public enum ColorOption {
        INSERT_SIZE,
        READ_STRAND,
        FIRST_OF_PAIR_STRAND,
        PAIR_ORIENTATION,
        SAMPLE,
        READ_GROUP,
        LIBRARY,
        MOVIE,
        ZMW,
        BISULFITE,
        NOMESEQ,
        TAG,
        NONE,
        UNEXPECTED_PAIR,
        MAPPED_SIZE,
        LINK_STRAND,
        YC_TAG,
        BASE_MODIFICATION
    }

    public enum SortOption {
        START, STRAND, NUCLEOTIDE, QUALITY, SAMPLE, READ_GROUP, INSERT_SIZE, FIRST_OF_PAIR_STRAND, MATE_CHR, TAG,
        SUPPLEMENTARY, NONE, HAPLOTYPE, READ_ORDER, READ_NAME
    }

    public enum GroupOption {
        STRAND("read strand"),
        SAMPLE("sample"),
        READ_GROUP("read group"),
        LIBRARY("library"),
        FIRST_OF_PAIR_STRAND("first-in-pair strand"),
        TAG("tag"),
        PAIR_ORIENTATION("pair orientation"),
        MATE_CHROMOSOME("chromosome of mate"),
        NONE("none"),
        SUPPLEMENTARY("supplementary flag"),
        BASE_AT_POS("base at position"),
        MOVIE("movie"),
        ZMW("ZMW"),
        HAPLOTYPE("haplotype"),
        READ_ORDER("read order"),
        LINKED("linked"),
        PHASE("phase"),
        REFERENCE_CONCORDANCE("reference concordance");

        public String label;

        GroupOption(String label) {
            this.label = label;
        }

    }

    public enum BisulfiteContext {
        CG, CHH, CHG, HCG, GCH, WCG, NONE
    }

    enum OrientationType {
        RR, LL, RL, LR, UNKNOWN
    }

    private static Logger log = LogManager.getLogger(AlignmentTrack.class);
    private static final int GROUP_LABEL_HEIGHT = 10;
    private static final int GROUP_MARGIN = 5;
    private static final int TOP_MARGIN = 20;
    private static final int DS_MARGIN_0 = 2;
    private static final int DOWNAMPLED_ROW_HEIGHT = 3;
    private static final int INSERTION_ROW_HEIGHT = 9;
    private static final int DS_MARGIN_2 = 5;
    private static int nClusters = 2;

    private static final Map<BisulfiteContext, String> bisulfiteContextToPubString = new HashMap<>();

    static {
        bisulfiteContextToPubString.put(BisulfiteContext.CG, "CG");
        bisulfiteContextToPubString.put(BisulfiteContext.CHH, "CHH");
        bisulfiteContextToPubString.put(BisulfiteContext.CHG, "CHG");
        bisulfiteContextToPubString.put(BisulfiteContext.HCG, "HCG");
        bisulfiteContextToPubString.put(BisulfiteContext.GCH, "GCH");
        bisulfiteContextToPubString.put(BisulfiteContext.WCG, "WCG");
        bisulfiteContextToPubString.put(BisulfiteContext.NONE, "None");
    }

    private static final Map<BisulfiteContext, Pair<byte[], byte[]>> bisulfiteContextToContextString = new HashMap<>();

    static {
        bisulfiteContextToContextString.put(BisulfiteContext.CG, new Pair<>(new byte[]{}, new byte[]{'G'}));
        bisulfiteContextToContextString.put(BisulfiteContext.CHH, new Pair<>(new byte[]{}, new byte[]{'H', 'H'}));
        bisulfiteContextToContextString.put(BisulfiteContext.CHG, new Pair<>(new byte[]{}, new byte[]{'H', 'G'}));
        bisulfiteContextToContextString.put(BisulfiteContext.HCG, new Pair<>(new byte[]{'H'}, new byte[]{'G'}));
        bisulfiteContextToContextString.put(BisulfiteContext.GCH, new Pair<>(new byte[]{'G'}, new byte[]{'H'}));
        bisulfiteContextToContextString.put(BisulfiteContext.WCG, new Pair<>(new byte[]{'W'}, new byte[]{'G'}));
    }

    public static boolean isBisulfiteColorType(ColorOption o) {
        return (o.equals(ColorOption.BISULFITE) || o.equals(ColorOption.NOMESEQ));
    }

    private static String getBisulfiteContextPubStr(BisulfiteContext item) {
        return bisulfiteContextToPubString.get(item);
    }

    public static byte[] getBisulfiteContextPreContext(BisulfiteContext item) {
        Pair<byte[], byte[]> pair = AlignmentTrack.bisulfiteContextToContextString.get(item);
        return pair.getFirst();
    }

    public static byte[] getBisulfiteContextPostContext(BisulfiteContext item) {
        Pair<byte[], byte[]> pair = AlignmentTrack.bisulfiteContextToContextString.get(item);
        return pair.getSecond();
    }


    private AlignmentDataManager dataManager;
    private SequenceTrack sequenceTrack;
    private CoverageTrack coverageTrack;
    private SpliceJunctionTrack spliceJunctionTrack;

    private final Genome genome;
    private ExperimentType experimentType;
    private final AlignmentRenderer renderer;
    RenderOptions renderOptions;

    private boolean removed = false;
    private RenderRollback renderRollback;
    private boolean showGroupLine;
    private Map<ReferenceFrame, List<InsertionInterval>> insertionIntervalsMap;
    private int expandedHeight = 14;
    private final int collapsedHeight = 9;
    private final int maxSquishedHeight = 5;
    private int squishedHeight = maxSquishedHeight;
    private final int minHeight = 50;

    private Rectangle alignmentsRect;
    private Rectangle downsampleRect;
    private Rectangle insertionRect;
    private ColorTable readNamePalette;

    // Dynamic fields
    protected final HashMap<String, Color> selectedReadNames = new HashMap<>();


    /**
     * Create a new alignment track
     *
     * @param locator
     * @param dataManager
     * @param genome
     */
    public AlignmentTrack(ResourceLocator locator, AlignmentDataManager dataManager, Genome genome) {
        super(locator);

        this.dataManager = dataManager;
        this.genome = genome;
        renderer = new AlignmentRenderer(this);
        renderOptions = new RenderOptions(this);
        setColor(DEFAULT_ALIGNMENT_COLOR);
        dataManager.setAlignmentTrack(this);
        dataManager.subscribe(this);

        IGVPreferences prefs = getPreferences();
        minimumHeight = 50;
        showGroupLine = prefs.getAsBoolean(SAM_SHOW_GROUP_SEPARATOR);
        try {
            setDisplayMode(DisplayMode.valueOf(prefs.get(SAM_DISPLAY_MODE).toUpperCase()));
        } catch (Exception e) {
            setDisplayMode(DisplayMode.EXPANDED);
        }
        if (prefs.getAsBoolean(SAM_SHOW_REF_SEQ)) {
            sequenceTrack = new SequenceTrack("Reference sequence");
            sequenceTrack.setHeight(14);
        }
        if (renderOptions.colorOption == ColorOption.BISULFITE) {
            setExperimentType(ExperimentType.BISULFITE);
        }
        readNamePalette = new PaletteColorTable(ColorUtilities.getDefaultPalette());
        insertionIntervalsMap = Collections.synchronizedMap(new HashMap<>());

        dataManager.setViewAsPairs(prefs.getAsBoolean(SAM_DISPLAY_PAIRED), renderOptions);

        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(AlignmentTrackEvent.class, this);
    }

    public void init() {
        if (experimentType == null) {
            ExperimentType type = dataManager.inferType();
            if (type != null) {
                setExperimentType(type);
            }
        }
    }


    @Override
    public void receiveEvent(Object event) {

        if (event instanceof FrameManager.ChangeEvent) {
            // Trim insertionInterval map to current frames
            Map<ReferenceFrame, List<InsertionInterval>> newMap = Collections.synchronizedMap(new HashMap<>());
            for (ReferenceFrame frame : ((FrameManager.ChangeEvent) event).getFrames()) {
                if (insertionIntervalsMap.containsKey(frame)) {
                    newMap.put(frame, insertionIntervalsMap.get(frame));
                }
            }
            insertionIntervalsMap = newMap;

        } else if (event instanceof AlignmentTrackEvent) {
            AlignmentTrackEvent e = (AlignmentTrackEvent) event;
            AlignmentTrackEvent.Type eventType = e.getType();
            switch (eventType) {
                case ALLELE_THRESHOLD:
                    dataManager.alleleThresholdChanged();
                    break;
                case RELOAD:
                    clearCaches();
                    repaint();
                case REFRESH:
                    repaint();
                    break;
            }

        }
    }

    void setExperimentType(ExperimentType type) {

        if (type != experimentType) {

            experimentType = type;

            boolean showJunction = getPreferences(type).getAsBoolean(Constants.SAM_SHOW_JUNCTION_TRACK);
            if (showJunction != spliceJunctionTrack.isVisible()) {
                spliceJunctionTrack.setVisible(showJunction);
                if (IGV.hasInstance()) {
                    IGV.getInstance().revalidateTrackPanels();
                }
            }

            boolean showCoverage = getPreferences(type).getAsBoolean(SAM_SHOW_COV_TRACK);
            if (showCoverage != coverageTrack.isVisible()) {
                coverageTrack.setVisible(showCoverage);
                if (IGV.hasInstance()) {
                    IGV.getInstance().revalidateTrackPanels();
                }
            }

            boolean showAlignments = getPreferences(type).getAsBoolean(SAM_SHOW_ALIGNMENT_TRACK);
            if (showAlignments != isVisible()) {
                setVisible(showAlignments);
                if (IGV.hasInstance()) {
                    IGV.getInstance().revalidateTrackPanels();
                }
            }
            //ExperimentTypeChangeEvent event = new ExperimentTypeChangeEvent(this, experimentType);
            //IGVEventBus.getInstance().post(event);
        }
    }

    ExperimentType getExperimentType() {
        return experimentType;
    }

    public AlignmentDataManager getDataManager() {
        return dataManager;
    }

    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public void setSpliceJunctionTrack(SpliceJunctionTrack spliceJunctionTrack) {
        this.spliceJunctionTrack = spliceJunctionTrack;
    }

    public SpliceJunctionTrack getSpliceJunctionTrack() {
        return spliceJunctionTrack;
    }

    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {
//
//        Alignment alignment = getAlignment(te);
//        if (alignment != null && alignment.getInsertions() != null) {
//            for (AlignmentBlock block : alignment.getInsertions()) {
//                if (block.containsPixel(te.getMouseEvent().getX())) {
//                    return new InsertionMenu(block);
//                }
//            }
//        }
        return new PopupMenu(te);
    }

    @Override
    public void setHeight(int preferredHeight) {
        super.setHeight(preferredHeight);
        minimumHeight = preferredHeight;
    }

    @Override
    public int getHeight() {

        int nGroups = dataManager.getMaxGroupCount();
        int h = Math.max(minHeight, getNLevels() * getRowHeight() + nGroups * GROUP_MARGIN + TOP_MARGIN
                + DS_MARGIN_0 + DOWNAMPLED_ROW_HEIGHT + DS_MARGIN_2);
        //if (insertionRect != null) {   // TODO - replace with expand insertions preference
        h += INSERTION_ROW_HEIGHT + DS_MARGIN_0;
        //}
        return Math.max(minimumHeight, h);
    }

    private int getRowHeight() {
        if (getDisplayMode() == DisplayMode.EXPANDED) {
            return expandedHeight;
        } else if (getDisplayMode() == DisplayMode.COLLAPSED) {
            return collapsedHeight;
        } else {
            return squishedHeight;
        }
    }

    private int getNLevels() {
        return dataManager.getNLevels();
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {

        if (frame.getChrName().equals(Globals.CHR_ALL) || frame.getScale() > dataManager.getMinVisibleScale()) {
            return true;   // Nothing to paint
        } else {
            List<InsertionInterval> insertionIntervals = getInsertionIntervals(frame);
            insertionIntervals.clear();
            return dataManager.isLoaded(frame);
        }
    }


    @Override
    public void load(ReferenceFrame referenceFrame) {
        if (log.isDebugEnabled()) {
            log.debug("Reading - thread: " + Thread.currentThread().getName());
        }
        dataManager.load(referenceFrame, renderOptions, true);
    }

    public void render(RenderContext context, Rectangle rect) {

        int viewWindowSize = context.getReferenceFrame().getCurrentRange().getLength();
        if (viewWindowSize > dataManager.getVisibilityWindow()) {
            Rectangle visibleRect = context.getVisibleRect().intersection(rect);
            Graphics2D g2 = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Zoom in to see alignments.", visibleRect, g2);
            return;
        }

        context.getGraphics2D("LABEL").setFont(FontManager.getFont(GROUP_LABEL_HEIGHT));

        // Split track rectangle into sections.
        int seqHeight = sequenceTrack == null ? 0 : sequenceTrack.getHeight();
        if (seqHeight > 0) {
            Rectangle seqRect = new Rectangle(rect);
            seqRect.height = seqHeight;
            sequenceTrack.render(context, seqRect);
        }

        // Top gap.
        rect.y += DS_MARGIN_0;

        downsampleRect = new Rectangle(rect);
        downsampleRect.height = DOWNAMPLED_ROW_HEIGHT;
        renderDownsampledIntervals(context, downsampleRect);

        if (renderOptions.isShowInsertionMarkers()) {
            insertionRect = new Rectangle(rect);
            insertionRect.y += DOWNAMPLED_ROW_HEIGHT + DS_MARGIN_0;
            insertionRect.height = INSERTION_ROW_HEIGHT;
            renderInsertionIntervals(context, insertionRect);
            rect.y = insertionRect.y + insertionRect.height;
        }

        alignmentsRect = new Rectangle(rect);
        alignmentsRect.y += 2;
        alignmentsRect.height -= (alignmentsRect.y - rect.y);
        renderAlignments(context, alignmentsRect);
    }

    private void renderDownsampledIntervals(RenderContext context, Rectangle downsampleRect) {

        // Might be offscreen
        if (!context.getVisibleRect().intersects(downsampleRect)) return;

        final AlignmentInterval loadedInterval = dataManager.getLoadedInterval(context.getReferenceFrame());
        if (loadedInterval == null) return;

        Graphics2D g = context.getGraphic2DForColor(Color.black);

        List<DownsampledInterval> intervals = loadedInterval.getDownsampledIntervals();
        for (DownsampledInterval interval : intervals) {
            final double scale = context.getScale();
            final double origin = context.getOrigin();

            int x0 = (int) ((interval.getStart() - origin) / scale);
            int x1 = (int) ((interval.getEnd() - origin) / scale);
            int w = Math.max(1, x1 - x0);
            // If there is room, leave a gap on one side
            if (w > 5) w--;
            // Greyscale from 0 -> 100 downsampled
            //int gray = 200 - interval.getCount();
            //Color color = (gray <= 0 ? Color.black : ColorUtilities.getGrayscaleColor(gray));
            g.fillRect(x0, downsampleRect.y, w, downsampleRect.height);
        }
    }


    private void renderAlignments(RenderContext context, Rectangle inputRect) {

        final AlignmentInterval loadedInterval = dataManager.getLoadedInterval(context.getReferenceFrame(), true);
        if (loadedInterval == null) {
            return;
        }

        final AlignmentCounts alignmentCounts = loadedInterval.getCounts();

        //log.debug("Render features");
        PackedAlignments groups = dataManager.getGroups(loadedInterval, renderOptions);
        if (groups == null) {
            //Assume we are still loading.
            //This might not always be true
            return;
        }

        // Check for YC tag
        if (renderOptions.colorOption == null && dataManager.hasYCTags()) {
            renderOptions.colorOption = ColorOption.YC_TAG;
        }

        Map<String, PEStats> peStats = dataManager.getPEStats();
        if (peStats != null) {
            renderOptions.peStats = peStats;
        }

        Rectangle visibleRect = context.getVisibleRect();

        // Divide rectangle into equal height levels
        double y = inputRect.getY();
        double h;
        if (getDisplayMode() == DisplayMode.EXPANDED) {
            h = expandedHeight;
        } else if (getDisplayMode() == DisplayMode.COLLAPSED) {
            h = collapsedHeight;
        } else {

            int visHeight = visibleRect.height;
            int depth = dataManager.getNLevels();
            if (depth == 0) {
                squishedHeight = Math.min(maxSquishedHeight, Math.max(1, expandedHeight));
            } else {
                squishedHeight = Math.min(maxSquishedHeight, Math.max(1, Math.min(expandedHeight, visHeight / depth)));
            }
            h = squishedHeight;
        }


        // Loop through groups
        Graphics2D groupBorderGraphics = context.getGraphic2DForColor(AlignmentRenderer.GROUP_DIVIDER_COLOR);
        int nGroups = groups.size();
        int groupNumber = 0;
        GroupOption groupOption = renderOptions.getGroupByOption();
        for (Map.Entry<String, List<Row>> entry : groups.entrySet()) {

            groupNumber++;
            double yGroup = y;  // Remember this for label

            // Loop through the alignment rows for this group
            List<Row> rows = entry.getValue();
            for (Row row : rows) {
                if ((visibleRect != null && y > visibleRect.getMaxY())) {
                    break;
                }

                assert visibleRect != null;
                if (y + h > visibleRect.getY()) {
                    Rectangle rowRectangle = new Rectangle(inputRect.x, (int) y, inputRect.width, (int) h);
                    renderer.renderAlignments(row.alignments, alignmentCounts, context, rowRectangle, renderOptions);
                    row.y = y;
                    row.h = h;
                }
                y += h;
            }

            if (groupOption != GroupOption.NONE) {
                // Draw a subtle divider line between groups
                if (showGroupLine) {
                    if (groupNumber < nGroups) {
                        int borderY = (int) y + GROUP_MARGIN / 2;
                        GraphicUtils.drawDottedDashLine(groupBorderGraphics, inputRect.x, borderY, inputRect.width, borderY);
                    }
                }

                // Label the group, if there is room
                double groupHeight = rows.size() * h;
                if (groupHeight > GROUP_LABEL_HEIGHT + 2) {
                    String groupName = entry.getKey();
                    Graphics2D g = context.getGraphics2D("LABEL");
                    FontMetrics fm = g.getFontMetrics();
                    Rectangle2D stringBouds = fm.getStringBounds(groupName, g);
                    Rectangle rect = new Rectangle(inputRect.x, (int) yGroup,
                            (int) stringBouds.getWidth() + 10, (int) stringBouds.getHeight());
                    GraphicUtils.drawVerticallyCenteredText(groupName, 5, rect, g, false, true);
                }
            }
            y += GROUP_MARGIN;
        }

        final int bottom = inputRect.y + inputRect.height;
        groupBorderGraphics.drawLine(inputRect.x, bottom, inputRect.width, bottom);
    }

    private List<InsertionInterval> getInsertionIntervals(ReferenceFrame frame) {
        List<InsertionInterval> insertionIntervals = insertionIntervalsMap.computeIfAbsent(frame, k -> new ArrayList<>());
        return insertionIntervals;
    }

    private void renderInsertionIntervals(RenderContext context, Rectangle rect) {

        // Might be offscreen
        if (!context.getVisibleRect().intersects(rect)) return;

        List<InsertionMarker> intervals = context.getInsertionMarkers();
        if (intervals == null) return;

        InsertionMarker selected = InsertionManager.getInstance().getSelectedInsertion(context.getChr());

        int w = (int) ((1.41 * rect.height) / 2);


        boolean hideSmallIndels = renderOptions.isHideSmallIndels();
        int smallIndelThreshold = renderOptions.getSmallIndelThreshold();

        List<InsertionInterval> insertionIntervals = getInsertionIntervals(context.getReferenceFrame());
        insertionIntervals.clear();
        for (InsertionMarker insertionMarker : intervals) {
            if (hideSmallIndels && insertionMarker.size < smallIndelThreshold) continue;

            final double scale = context.getScale();
            final double origin = context.getOrigin();
            int midpoint = (int) ((insertionMarker.position - origin) / scale);
            int x0 = midpoint - w;
            int x1 = midpoint + w;

            Rectangle iRect = new Rectangle(x0 + context.translateX, rect.y, 2 * w, rect.height);

            insertionIntervals.add(new InsertionInterval(iRect, insertionMarker));

            Color c = (selected != null && selected.position == insertionMarker.position) ? new Color(200, 0, 0, 80) : AlignmentRenderer.purple;
            Graphics2D g = context.getGraphic2DForColor(c);


            g.fillPolygon(new Polygon(new int[]{x0, x1, midpoint},
                    new int[]{rect.y, rect.y, rect.y + rect.height}, 3));
        }
    }

    public void renderExpandedInsertion(InsertionMarker insertionMarker, RenderContext context, Rectangle inputRect) {


        boolean leaveMargin = getDisplayMode() != DisplayMode.SQUISHED;

        // Insertion interval
        Graphics2D g = context.getGraphic2DForColor(Color.red);
        Rectangle iRect = new Rectangle(inputRect.x, insertionRect.y, inputRect.width, insertionRect.height);
        g.fill(iRect);
        List<InsertionInterval> insertionIntervals = getInsertionIntervals(context.getReferenceFrame());

        iRect.x += context.translateX;
        insertionIntervals.add(new InsertionInterval(iRect, insertionMarker));


        inputRect.y += DS_MARGIN_0 + DOWNAMPLED_ROW_HEIGHT + DS_MARGIN_0 + INSERTION_ROW_HEIGHT + DS_MARGIN_2;

        //log.debug("Render features");
        final AlignmentInterval loadedInterval = dataManager.getLoadedInterval(context.getReferenceFrame(), true);
        PackedAlignments groups = dataManager.getGroups(loadedInterval, renderOptions);
        if (groups == null) {
            //Assume we are still loading.
            //This might not always be true
            return;
        }

        Rectangle visibleRect = context.getVisibleRect();

        // Divide rectangle into equal height levels
        double y = inputRect.getY() - 3;
        double h;
        if (getDisplayMode() == DisplayMode.EXPANDED) {
            h = expandedHeight;
        } else if (getDisplayMode() == DisplayMode.COLLAPSED) {
            h = collapsedHeight;
        } else {
            int visHeight = visibleRect.height;
            int depth = dataManager.getNLevels();
            if (depth == 0) {
                squishedHeight = Math.min(maxSquishedHeight, Math.max(1, expandedHeight));
            } else {
                squishedHeight = Math.min(maxSquishedHeight, Math.max(1, Math.min(expandedHeight, visHeight / depth)));
            }
            h = squishedHeight;
        }


        for (Map.Entry<String, List<Row>> entry : groups.entrySet()) {


            // Loop through the alignment rows for this group
            List<Row> rows = entry.getValue();
            for (Row row : rows) {
                if ((visibleRect != null && y > visibleRect.getMaxY())) {
                    return;
                }

                assert visibleRect != null;
                if (y + h > visibleRect.getY()) {
                    Rectangle rowRectangle = new Rectangle(inputRect.x, (int) y, inputRect.width, (int) h);
                    renderer.renderExpandedInsertion(insertionMarker, row.alignments, context, rowRectangle, leaveMargin);
                    row.y = y;
                    row.h = h;
                }
                y += h;
            }

            y += GROUP_MARGIN;


        }

    }

    private InsertionInterval getInsertionInterval(ReferenceFrame frame, int x, int y) {
        List<InsertionInterval> insertionIntervals = getInsertionIntervals(frame);
        for (InsertionInterval i : insertionIntervals) {
            if (i.rect.contains(x, y)) return i;
        }

        return null;
    }

    /**
     * Sort alignment rows based on alignments that intersect location
     *
     * @return Whether sorting was performed. If data is still loading, this will return false
     */
    public boolean sortRows(SortOption option, ReferenceFrame referenceFrame, double location, String tag) {
        return dataManager.sortRows(option, referenceFrame, location, tag);
    }

    private static void sortAlignmentTracks(SortOption option, String tag) {
        IGV.getInstance().sortAlignmentTracks(option, tag);
        Collection<IGVPreferences> allPrefs = PreferencesManager.getAllPreferences();
        for (IGVPreferences prefs : allPrefs) {
            prefs.put(SAM_SORT_OPTION, option.toString());
            prefs.put(SAM_SORT_BY_TAG, tag);
        }
    }

    /**
     * Visually regroup alignments by the provided {@code GroupOption}.
     *
     * @param option
     * @see AlignmentDataManager#packAlignments
     */
    public void groupAlignments(GroupOption option, String tag, Range pos) {
        if (option == GroupOption.TAG && tag != null) {
            renderOptions.setGroupByTag(tag);
        }
        if (option == GroupOption.BASE_AT_POS && pos != null) {
            renderOptions.setGroupByPos(pos);
        }
        renderOptions.setGroupByOption(option);
        dataManager.packAlignments(renderOptions);
    }

    public void setBisulfiteContext(BisulfiteContext option) {
        renderOptions.bisulfiteContext = option;
        getPreferences().put(SAM_BISULFITE_CONTEXT, option.toString());
    }

    public void setColorOption(ColorOption option) {
        renderOptions.setColorOption(option);
    }

    public void setColorByTag(String tag) {
        renderOptions.setColorByTag(tag);
        getPreferences(experimentType).put(SAM_COLOR_BY_TAG, tag);
    }

    public void packAlignments() {
        dataManager.packAlignments(renderOptions);
    }

    /**
     * Copy the contents of the popup text to the system clipboard.
     */
    private void copyToClipboard(final TrackClickEvent e, Alignment alignment, double location, int mouseX) {

        if (alignment != null) {
            StringBuilder buf = new StringBuilder();
            buf.append(alignment.getClipboardString(location, mouseX)
                    .replace("<b>", "")
                    .replace("</b>", "")
                    .replace("<br>", "\n")
                    .replace("<br/>", "\n")
                    .replace("<hr>", "\n------------------\n")
                    .replace("<hr/>", "\n------------------\n"));
            buf.append("\n");
            buf.append("Alignment start position = ").append(alignment.getChr()).append(":").append(alignment.getAlignmentStart() + 1);
            buf.append("\n");
            buf.append(alignment.getReadSequence());
            StringSelection stringSelection = new StringSelection(buf.toString());
            Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
            clipboard.setContents(stringSelection, null);
        }

    }

    /**
     * Jump to the mate region
     */
    private void gotoMate(final TrackClickEvent te, Alignment alignment) {


        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {

                setSelected(alignment);

                String chr = mate.getChr();
                int start = mate.start - 1;

                // Don't change scale
                double range = te.getFrame().getEnd() - te.getFrame().getOrigin();
                int newStart = (int) Math.max(0, (start + (alignment.getEnd() - alignment.getStart()) / 2 - range / 2));
                int newEnd = newStart + (int) range;
                te.getFrame().jumpTo(chr, newStart, newEnd);
                te.getFrame().recordHistory();
            } else {
                MessageUtils.showMessage("Alignment does not have mate, or it is not mapped.");
            }
        }
    }

    /**
     * Split the screen so the current view and mate region are side by side.
     * Need a better name for this method.
     */
    private void splitScreenMate(final TrackClickEvent te, Alignment alignment) {

        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {

                setSelected(alignment);

                String mateChr = mate.getChr();
                int mateStart = mate.start - 1;

                ReferenceFrame frame = te.getFrame();
                String locus1 = frame.getFormattedLocusString();

                // Generate a locus string for the read mate.  Keep the window width (in base pairs) == to the current range
                Range range = frame.getCurrentRange();
                int length = range.getLength();
                int s2 = Math.max(0, mateStart - length / 2);
                int e2 = s2 + length;
                String startStr = String.valueOf(s2);
                String endStr = String.valueOf(e2);
                String mateLocus = mateChr + ":" + startStr + "-" + endStr;

                Session currentSession = IGV.getInstance().getSession();

                List<String> loci;
                if (FrameManager.isGeneListMode()) {
                    loci = new ArrayList<>(FrameManager.getFrames().size());
                    for (ReferenceFrame ref : FrameManager.getFrames()) {
                        //If the frame-name is a locus, we use it unaltered
                        //Don't want to reprocess, easy to get off-by-one
                        String name = ref.getName();
                        if (Locus.fromString(name) != null) {
                            loci.add(name);
                        } else {
                            loci.add(ref.getFormattedLocusString());
                        }

                    }
                    loci.add(mateLocus);
                } else {
                    loci = Arrays.asList(locus1, mateLocus);
                }

                StringBuilder listName = new StringBuilder();
                for (String s : loci) {
                    listName.append(s + "   ");
                }

                GeneList geneList = new GeneList(listName.toString(), loci, false);
                currentSession.setCurrentGeneList(geneList);

                Comparator<String> geneListComparator = (n0, n1) -> {
                    ReferenceFrame f0 = FrameManager.getFrame(n0);
                    ReferenceFrame f1 = FrameManager.getFrame(n1);

                    String chr0 = f0 == null ? "" : f0.getChrName();
                    String chr1 = f1 == null ? "" : f1.getChrName();
                    int s0 = f0 == null ? 0 : f0.getCurrentRange().getStart();
                    int s1 = f1 == null ? 0 : f1.getCurrentRange().getStart();

                    int chrComp = ChromosomeNameComparator.get().compare(chr0, chr1);
                    if (chrComp != 0) return chrComp;
                    return s0 - s1;
                };

                //Need to sort the frames by position
                currentSession.sortGeneList(geneListComparator);
                IGV.getInstance().resetFrames();
            } else {
                MessageUtils.showMessage("Alignment does not have mate, or it is not mapped.");
            }
        }
    }

    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0.0f;
    }


    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        if (downsampleRect != null && mouseY > downsampleRect.y && mouseY <= downsampleRect.y + downsampleRect.height) {
            AlignmentInterval loadedInterval = dataManager.getLoadedInterval(frame);
            if (loadedInterval == null) {
                return null;
            } else {
                List<DownsampledInterval> intervals = loadedInterval.getDownsampledIntervals();
                DownsampledInterval interval = FeatureUtils.getFeatureAt(position, 0, intervals);
                if (interval != null) {
                    return interval.getValueString();
                }
                return null;
            }
        } else {

            InsertionInterval insertionInterval = getInsertionInterval(frame, mouseX, mouseY);
            if (insertionInterval != null) {
                return "Insertions (" + insertionInterval.insertionMarker.size + " bases)";
            } else {
                Alignment feature = getAlignmentAt(position, mouseY, frame);
                if (feature != null) {
                    return feature.getAlignmentValueString(position, mouseX, renderOptions);
                }
            }

        }
        return null;
    }


    private Alignment getAlignment(final TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();
        final ReferenceFrame frame = te.getFrame();
        if (frame == null) {
            return null;
        }
        final double location = frame.getChromosomePosition(e.getX());
        return getAlignmentAt(location, e.getY(), frame);
    }

    private Alignment getAlignmentAt(double position, int y, ReferenceFrame frame) {

        if (alignmentsRect == null || dataManager == null) {
            return null;   // <= not loaded yet
        }
        PackedAlignments groups = dataManager.getGroupedAlignmentsContaining(position, frame);

        if (groups == null || groups.isEmpty()) {
            return null;
        }

        for (List<Row> rows : groups.values()) {
            for (Row row : rows) {
                if (y >= row.y && y <= row.y + row.h) {
                    List<Alignment> features = row.alignments;

                    // No buffer for alignments,  you must zoom in far enough for them to be visible
                    int buffer = 0;
                    return FeatureUtils.getFeatureAt(position, buffer, features);
                }
            }
        }
        return null;
    }


    /**
     * Get the most "specific" alignment at the specified location.  Specificity refers to the smallest alignemnt
     * in a group that contains the location (i.e. if a group of linked alignments overlap take the smallest one).
     *
     * @param te
     * @return
     */
    private Alignment getSpecficAlignment(TrackClickEvent te) {

        Alignment alignment = getAlignment(te);
        if (alignment != null) {
            final ReferenceFrame frame = te.getFrame();
            MouseEvent e = te.getMouseEvent();
            final double location = frame.getChromosomePosition(e.getX());

            if (alignment instanceof LinkedAlignment) {

                Alignment sa = null;
                for (Alignment a : ((LinkedAlignment) alignment).alignments) {
                    if (a.contains(location)) {
                        if (sa == null || (a.getAlignmentEnd() - a.getAlignmentStart() < sa.getAlignmentEnd() - sa.getAlignmentStart())) {
                            sa = a;
                        }
                    }
                }
                alignment = sa;

            } else if (alignment instanceof PairedAlignment) {
                Alignment sa = null;
                if (((PairedAlignment) alignment).firstAlignment.contains(location)) {
                    sa = ((PairedAlignment) alignment).firstAlignment;
                } else if (((PairedAlignment) alignment).secondAlignment.contains(location)) {
                    sa = ((PairedAlignment) alignment).secondAlignment;
                }
                alignment = sa;
            }
        }
        return alignment;
    }


    @Override
    public boolean handleDataClick(TrackClickEvent te) {

        MouseEvent e = te.getMouseEvent();
        if (Globals.IS_MAC && e.isMetaDown() || (!Globals.IS_MAC && e.isControlDown())) {
            // Selection
            final ReferenceFrame frame = te.getFrame();
            if (frame != null) {
                selectAlignment(e, frame);
                IGV.getInstance().repaint(this);
                return true;
            }
        }

        InsertionInterval insertionInterval = getInsertionInterval(te.getFrame(), te.getMouseEvent().getX(), te.getMouseEvent().getY());
        if (insertionInterval != null) {

            final String chrName = te.getFrame().getChrName();
            InsertionMarker currentSelection = InsertionManager.getInstance().getSelectedInsertion(chrName);
            if (currentSelection != null && currentSelection.position == insertionInterval.insertionMarker.position) {
                InsertionManager.getInstance().clearSelected();
            } else {
                InsertionManager.getInstance().setSelected(chrName, insertionInterval.insertionMarker.position);
            }

            IGVEventBus.getInstance().post(new InsertionSelectionEvent(insertionInterval.insertionMarker));

            return true;
        }


        if (IGV.getInstance().isShowDetailsOnClick()) {
            openTooltipWindow(te);
            return true;
        }

        return false;
    }

    private void selectAlignment(MouseEvent e, ReferenceFrame frame) {
        double location = frame.getChromosomePosition(e.getX());
        Alignment alignment = this.getAlignmentAt(location, e.getY(), frame);
        if (alignment != null) {
            if (selectedReadNames.containsKey(alignment.getReadName())) {
                selectedReadNames.remove(alignment.getReadName());
            } else {
                setSelected(alignment);
            }

        }
    }

    private void setSelected(Alignment alignment) {
        Color c = readNamePalette.get(alignment.getReadName());
        selectedReadNames.put(alignment.getReadName(), c);
    }

    private void clearCaches() {
        if (dataManager != null) dataManager.clear();
        if (spliceJunctionTrack != null) spliceJunctionTrack.clear();
    }


    public void setViewAsPairs(boolean vAP) {
        // TODO -- generalize this test to all incompatible pairings
        if (vAP && renderOptions.groupByOption == GroupOption.STRAND) {
            boolean ungroup = MessageUtils.confirm("\"View as pairs\" is incompatible with \"Group by strand\". Ungroup?");
            if (ungroup) {
                renderOptions.setGroupByOption(null);
            } else {
                return;
            }
        }

        dataManager.setViewAsPairs(vAP, renderOptions);
        repaint();
    }


    public enum ExperimentType {OTHER, RNA, BISULFITE, THIRD_GEN}


    class RenderRollback {
        final ColorOption colorOption;
        final GroupOption groupByOption;
        final String groupByTag;
        final String colorByTag;
        final String linkByTag;
        final DisplayMode displayMode;
        final int expandedHeight;
        final boolean showGroupLine;

        RenderRollback(RenderOptions renderOptions, DisplayMode displayMode) {
            this.colorOption = renderOptions.colorOption;
            this.groupByOption = renderOptions.groupByOption;
            this.colorByTag = renderOptions.colorByTag;
            this.groupByTag = renderOptions.groupByTag;
            this.displayMode = displayMode;
            this.expandedHeight = AlignmentTrack.this.expandedHeight;
            this.showGroupLine = AlignmentTrack.this.showGroupLine;
            this.linkByTag = renderOptions.linkByTag;
        }

        void restore(RenderOptions renderOptions) {
            renderOptions.colorOption = this.colorOption;
            renderOptions.groupByOption = this.groupByOption;
            renderOptions.colorByTag = this.colorByTag;
            renderOptions.groupByTag = this.groupByTag;
            renderOptions.linkByTag = this.linkByTag;
            AlignmentTrack.this.expandedHeight = this.expandedHeight;
            AlignmentTrack.this.showGroupLine = this.showGroupLine;
            AlignmentTrack.this.setDisplayMode(this.displayMode);
        }
    }

    public boolean isRemoved() {
        return removed;
    }

    @Override
    public boolean isVisible() {
        return super.isVisible() && !removed;
    }

    IGVPreferences getPreferences() {
        return getPreferences(experimentType);
    }

    public static IGVPreferences getPreferences(ExperimentType type) {

        try {
            // Disable experimentType preferences for 2.4
            if (Globals.VERSION.contains("2.4")) {
                return PreferencesManager.getPreferences(NULL_CATEGORY);
            } else {
                String prefKey = Constants.NULL_CATEGORY;
                if (type == ExperimentType.THIRD_GEN) {
                    prefKey = Constants.THIRD_GEN;
                } else if (type == ExperimentType.RNA) {
                    prefKey = Constants.RNA;
                }
                return PreferencesManager.getPreferences(prefKey);
            }
        } catch (NullPointerException e) {
            String prefKey = Constants.NULL_CATEGORY;
            if (type == ExperimentType.THIRD_GEN) {
                prefKey = Constants.THIRD_GEN;
            } else if (type == ExperimentType.RNA) {
                prefKey = Constants.RNA;
            }
            return PreferencesManager.getPreferences(prefKey);
        }
    }


    @Override
    public void unload() {
        super.unload();
        if (dataManager != null) {
            dataManager.unsubscribe(this);
        }
        removed = true;
        setVisible(false);
    }

    private boolean isLinkedReads() {
        return renderOptions != null && renderOptions.isLinkedReads();
    }

    /**
     * Set the view to a 10X style "linked-read view".  This option tries to achieve a view similar to the
     * 10X Loupe view described here: https://support.10xgenomics.com/genome-exome/software/visualization/latest/linked-reads
     *
     * @param linkedReads
     * @param tag
     */
    private void setLinkedReadView(boolean linkedReads, String tag) {
        if (!linkedReads || isLinkedReadView()) {
            undoLinkedReadView();
        }
        renderOptions.setLinkedReads(linkedReads);
        if (linkedReads) {
            renderOptions.setLinkByTag(tag);
            renderOptions.setColorOption(ColorOption.TAG);
            renderOptions.setColorByTag(tag);
            if (dataManager.isPhased()) {
                renderOptions.setGroupByOption(GroupOption.TAG);
                renderOptions.setGroupByTag("HP");
            }
            showGroupLine = false;
            setDisplayMode(DisplayMode.SQUISHED);
        }
        dataManager.packAlignments(renderOptions);
        repaint();
    }

    /**
     * Detect if we are in linked-read view
     */
    private boolean isLinkedReadView() {
        return renderOptions != null &&
                renderOptions.isLinkedReads() &&
                renderOptions.getLinkByTag() != null &&
                renderOptions.getColorOption() == ColorOption.TAG &&
                renderOptions.getColorByTag() != null;
    }

    /**
     * Link alignments by arbitrary tag, without the extra settings applied to link-read-view
     *
     * @param linkReads
     * @param tag
     */
    private void setLinkByTag(boolean linkReads, String tag) {
        if (isLinkedReadView()) {
            undoLinkedReadView();
        }
        if (linkReads) {
            renderOptions.setLinkByTag(tag);
            if (renderOptions.getGroupByOption() == GroupOption.NONE) {
                renderOptions.setGroupByOption(GroupOption.LINKED);
            }
        } else {
            renderOptions.setLinkByTag(null);
            if (renderOptions.getGroupByOption() == GroupOption.LINKED) {
                renderOptions.setGroupByOption(GroupOption.NONE);
            }
        }
        renderOptions.setLinkedReads(linkReads);
        dataManager.packAlignments(renderOptions);
        repaint();
    }

    private void undoLinkedReadView() {
        renderOptions.setLinkByTag(null);
        renderOptions.setColorOption(ColorOption.NONE);
        renderOptions.setColorByTag(null);
        renderOptions.setGroupByOption(GroupOption.NONE);
        renderOptions.setGroupByTag(null);
        showGroupLine = true;
        setDisplayMode(DisplayMode.EXPANDED);
    }

    private void sendPairsToCircularView(TrackClickEvent e) {

        List<ReferenceFrame> frames = e.getFrame() != null ?
                Arrays.asList(e.getFrame()) :
                FrameManager.getFrames();

        List<Alignment> inView = new ArrayList<>();
        for (ReferenceFrame frame : frames) {
            AlignmentInterval interval = AlignmentTrack.this.getDataManager().getLoadedInterval(frame);
            if (interval != null) {
                Iterator<Alignment> iter = interval.getAlignmentIterator();
                Range r = frame.getCurrentRange();
                while (iter.hasNext()) {
                    Alignment a = iter.next();
                    if (a.getEnd() > r.getStart() && a.getStart() < r.getEnd()) {
                        final boolean isDiscordantPair = a.isPaired() && a.getMate().isMapped() &&
                                (!a.getMate().getChr().equals(a.getChr()) || Math.abs(a.getInferredInsertSize()) > 10000);
                        if (isDiscordantPair) {
                            inView.add(a);
                        }
                    }
                }
            }
            Color chordColor = AlignmentTrack.this.getColor().equals(DEFAULT_ALIGNMENT_COLOR) ? Color.BLUE : AlignmentTrack.this.getColor();
            CircularViewUtilities.sendAlignmentsToJBrowse(inView, AlignmentTrack.this.getName(), chordColor);
        }
    }

    private void sendSplitToCircularView(TrackClickEvent e) {

        List<ReferenceFrame> frames = e.getFrame() != null ?
                Arrays.asList(e.getFrame()) :
                FrameManager.getFrames();

        List<Alignment> inView = new ArrayList<>();
        for (ReferenceFrame frame : frames) {
            AlignmentInterval interval = AlignmentTrack.this.getDataManager().getLoadedInterval(frame);
            if (interval != null) {
                Iterator<Alignment> iter = interval.getAlignmentIterator();
                Range r = frame.getCurrentRange();
                while (iter.hasNext()) {
                    Alignment a = iter.next();
                    if (a.getEnd() > r.getStart() && a.getStart() < r.getEnd() && a.getAttribute("SA") != null) {
                        inView.add(a);
                    }
                }
            }
            Color chordColor = AlignmentTrack.this.getColor().equals(DEFAULT_ALIGNMENT_COLOR) ? Color.BLUE : AlignmentTrack.this.getColor();
            CircularViewUtilities.sendAlignmentsToJBrowse(inView, AlignmentTrack.this.getName(), chordColor);
        }
    }

    /**
     * Listener for deselecting one component when another is selected
     */
    private static class Deselector implements ActionListener {

        private final JMenuItem toDeselect;
        private final JMenuItem parent;

        Deselector(JMenuItem parent, JMenuItem toDeselect) {
            this.parent = parent;
            this.toDeselect = toDeselect;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            if (this.parent.isSelected()) {
                this.toDeselect.setSelected(false);
            }
        }
    }

    private static class InsertionInterval {

        final Rectangle rect;
        final InsertionMarker insertionMarker;

        InsertionInterval(Rectangle rect, InsertionMarker insertionMarker) {
            this.rect = rect;
            this.insertionMarker = insertionMarker;
        }
    }

    class PopupMenu extends IGVPopupMenu {


        PopupMenu(final TrackClickEvent e) {

            final MouseEvent me = e.getMouseEvent();
            ReferenceFrame frame = e.getFrame();
            Alignment clickedAlignment = null;

            if (frame != null) {
                double location = frame.getChromosomePosition(me.getX());
                clickedAlignment = getAlignmentAt(location, me.getY(), frame);
            }


            Collection<Track> tracks = new ArrayList();
            tracks.add(AlignmentTrack.this);

            JLabel popupTitle = new JLabel("  " + AlignmentTrack.this.getName(), JLabel.CENTER);

            Font newFont = getFont().deriveFont(Font.BOLD, 12);
            popupTitle.setFont(newFont);
            add(popupTitle);

            if (PreferencesManager.getPreferences().getAsBoolean(CIRC_VIEW_ENABLED) && CircularViewUtilities.ping()) {
                addSeparator();
                JMenuItem item = new JMenuItem("Add Discordant Pairs to Circular View");
                //item.setEnabled(dataManager.isPairedEnd());
                add(item);
                item.addActionListener(ae -> AlignmentTrack.this.sendPairsToCircularView(e));

                JMenuItem item2 = new JMenuItem("Add Split Reads to Circular View");
                //item.setEnabled(dataManager.isPairedEnd());
                add(item2);
                item2.addActionListener(ae -> AlignmentTrack.this.sendSplitToCircularView(e));
            }

            addSeparator();
            add(TrackMenuUtils.getTrackRenameItem(tracks));
            addCopyToClipboardItem(e, clickedAlignment);

            //         addSeparator();
            //          addExpandInsertions();

            addSeparator();
            JMenuItem item = new JMenuItem("Change Track Color...");
            item.addActionListener(evt -> TrackMenuUtils.changeTrackColor(tracks));
            add(item);

            addSeparator();
            addExperimentTypeMenuItem();

            if (experimentType == ExperimentType.THIRD_GEN) {
                addHaplotype(e);
            }

            addLinkedReadItems();

            addSeparator();
            addGroupMenuItem(e);
            addSortMenuItem();
            addColorByMenuItem();
            //addFilterMenuItem();
            addPackMenuItem();

            addSeparator();
            addShadeBaseByMenuItem();
            JMenuItem misMatchesItem = addShowMismatchesMenuItem();
            JMenuItem showAllItem = addShowAllBasesMenuItem();

            misMatchesItem.addActionListener(new Deselector(misMatchesItem, showAllItem));
            showAllItem.addActionListener(new Deselector(showAllItem, misMatchesItem));

            // Paired end items
            addSeparator();
            addViewAsPairsMenuItem();
            if (clickedAlignment != null) {
                addGoToMate(e, clickedAlignment);
                showMateRegion(e, clickedAlignment);
            }
            addInsertSizeMenuItem();

            // Third gen (primarily) items
            addSeparator();
            addThirdGenItems();

            // Display mode items
            addSeparator();
            TrackMenuUtils.addDisplayModeItems(tracks, this);

            // Select items
            addSeparator();
            addSelectByNameItem();
            addClearSelectionsMenuItem();

            // Copy items
            addSeparator();
            addCopySequenceItems(e);
            addConsensusSequence(e);

            // Blat items
            addSeparator();
            addBlatItem(e);
            addBlatClippingItems(e);

            AlignmentBlock insertion = getInsertion(clickedAlignment, e.getMouseEvent().getX());
            if (insertion != null) {
                addSeparator();
                addInsertionItems(insertion);
            }

            addSeparator();
            JMenuItem sashimi = new JMenuItem("Sashimi Plot");
            sashimi.addActionListener(e1 -> SashimiPlot.getSashimiPlot(null));
            add(sashimi);

            addSeparator();
            addShowItems();

//
//            if (getPreferences().get(Constants.EXTVIEW_URL) != null) {
//                addSeparator();
//                addExtViewItem(e);
//            }

        }


        private void addHaplotype(TrackClickEvent e) {

            JMenuItem item = new JMenuItem("Cluster (phase) alignments");

            final ReferenceFrame frame;
            if (e.getFrame() == null && FrameManager.getFrames().size() == 1) {
                frame = FrameManager.getFrames().get(0);
            } else {
                frame = e.getFrame();
            }

            item.setEnabled(frame != null);
            add(item);

            item.addActionListener(ae -> {
                //This shouldn't ever be true, but just in case it's more user-friendly
                if (frame == null) {
                    MessageUtils.showMessage("Unknown region bounds");
                    return;
                }

                String nString = MessageUtils.showInputDialog("Enter the number of clusters", String.valueOf(AlignmentTrack.nClusters));
                if (nString == null) {
                    return;
                }
                try {
                    AlignmentTrack.nClusters = Integer.parseInt(nString);
                } catch (NumberFormatException e1) {
                    MessageUtils.showMessage("Clusters size must be an integer");
                    return;
                }

                final int start = (int) frame.getOrigin();
                final int end = (int) frame.getEnd();

                AlignmentInterval interval = dataManager.getLoadedInterval(frame);
                HaplotypeUtils haplotypeUtils = new HaplotypeUtils(interval, AlignmentTrack.this.genome);
                boolean success = haplotypeUtils.clusterAlignments(frame.getChrName(), start, end, AlignmentTrack.nClusters);

                if (success) {
                    AlignmentTrack.this.groupAlignments(GroupOption.HAPLOTYPE, null, null);
                    AlignmentTrack.this.repaint();
                }

                //dataManager.sortRows(SortOption.HAPLOTYPE, frame, (end + start) / 2, null);
                //AlignmentTrack.repaint();

            });


        }


        public JMenuItem addExpandInsertions() {

            final JMenuItem item = new JCheckBoxMenuItem("Expand insertions");
            final Session session = IGV.getInstance().getSession();
            item.setSelected(session.expandInsertions);

            item.addActionListener(aEvt -> {
                session.expandInsertions = !session.expandInsertions;
                AlignmentTrack.this.repaint();
            });
            add(item);
            return item;
        }

        /**
         * Item for exporting "consensus" sequence of region, based
         * on loaded alignments.
         *
         * @param e
         */
        private void addConsensusSequence(TrackClickEvent e) {
            //Export consensus sequence
            JMenuItem item = new JMenuItem("Copy consensus sequence");


            final ReferenceFrame frame;
            if (e.getFrame() == null && FrameManager.getFrames().size() == 1) {
                frame = FrameManager.getFrames().get(0);
            } else {
                frame = e.getFrame();
            }

            item.setEnabled(frame != null);
            add(item);

            item.addActionListener(ae -> {
                //This shouldn't ever be true, but just in case it's more user-friendly
                if (frame == null) {
                    MessageUtils.showMessage("Unknown region bounds, cannot export consensus");
                    return;
                }
                final int start = (int) frame.getOrigin();
                final int end = (int) frame.getEnd();
                if ((end - start) > 1000000) {
                    MessageUtils.showMessage("Cannot export region more than 1 Megabase");
                    return;
                }
                AlignmentInterval interval = dataManager.getLoadedInterval(frame);
                AlignmentCounts counts = interval.getCounts();
                String text = PFMExporter.createPFMText(counts, frame.getChrName(), start, end);
                StringUtils.copyTextToClipboard(text);
            });


        }

        private JMenu getBisulfiteContextMenuItem(ButtonGroup group) {
            // Change track height by attribute
            //JMenu bisulfiteContextMenu = new JMenu("Bisulfite Contexts");
            JMenu bisulfiteContextMenu = new JMenu("bisulfite mode");


            JRadioButtonMenuItem nomeESeqOption = null;
            boolean showNomeESeq = getPreferences().getAsBoolean(SAM_NOMESEQ_ENABLED);
            if (showNomeESeq) {
                nomeESeqOption = new JRadioButtonMenuItem("NOMe-seq bisulfite mode");
                nomeESeqOption.setSelected(renderOptions.getColorOption() == ColorOption.NOMESEQ);
                nomeESeqOption.addActionListener(aEvt -> {
                    setColorOption(ColorOption.NOMESEQ);
                    AlignmentTrack.this.repaint();
                });
                group.add(nomeESeqOption);
            }

            for (final BisulfiteContext item : BisulfiteContext.values()) {

                String optionStr = getBisulfiteContextPubStr(item);
                JRadioButtonMenuItem m1 = new JRadioButtonMenuItem(optionStr);
                m1.setSelected(renderOptions.bisulfiteContext == item);
                m1.addActionListener(aEvt -> {
                    setColorOption(ColorOption.BISULFITE);
                    setBisulfiteContext(item);
                    AlignmentTrack.this.repaint();
                });
                bisulfiteContextMenu.add(m1);
                group.add(m1);
            }

            if (nomeESeqOption != null) {
                bisulfiteContextMenu.add(nomeESeqOption);
            }

            return bisulfiteContextMenu;

        }

        void addSelectByNameItem() {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("Select by name...");
            item.addActionListener(aEvt -> {
                String val = MessageUtils.showInputDialog("Enter read name: ");
                if (val != null && val.trim().length() > 0) {
                    selectedReadNames.put(val, readNamePalette.get(val));
                    AlignmentTrack.this.repaint();
                }
            });
            add(item);
        }

        void addExperimentTypeMenuItem() {
            Map<String, ExperimentType> mappings = new LinkedHashMap<>();
            mappings.put("Other", ExperimentType.OTHER);
            mappings.put("RNA", ExperimentType.RNA);
            mappings.put("3rd Gen", ExperimentType.THIRD_GEN);
            //mappings.put("Bisulfite", ExperimentType.BISULFITE);
            JMenu groupMenu = new JMenu("Experiment Type");
            ButtonGroup group = new ButtonGroup();
            for (Map.Entry<String, ExperimentType> el : mappings.entrySet()) {
                JCheckBoxMenuItem mi = getExperimentTypeMenuItem(el.getKey(), el.getValue());
                groupMenu.add(mi);
                group.add(mi);
            }
            add(groupMenu);
        }

        private JCheckBoxMenuItem getExperimentTypeMenuItem(String label, final ExperimentType option) {
            JCheckBoxMenuItem mi = new JCheckBoxMenuItem(label);
            mi.setSelected(AlignmentTrack.this.getExperimentType() == option);
            mi.addActionListener(aEvt -> AlignmentTrack.this.setExperimentType(option));
            return mi;
        }

        void addGroupMenuItem(final TrackClickEvent te) {//ReferenceFrame frame) {
            final MouseEvent me = te.getMouseEvent();
            ReferenceFrame frame = te.getFrame();
            if (frame == null) {
                frame = FrameManager.getDefaultFrame();  // Clicked over name panel, not a specific frame
            }
            final Range range = frame.getCurrentRange();
            final String chrom = range.getChr();
            final int chromStart = (int) frame.getChromosomePosition(me.getX());
            // Change track height by attribute
            JMenu groupMenu = new JMenu("Group alignments by");
            ButtonGroup group = new ButtonGroup();

            GroupOption[] groupOptions = {
                    GroupOption.NONE, GroupOption.STRAND, GroupOption.FIRST_OF_PAIR_STRAND, GroupOption.SAMPLE,
                    GroupOption.LIBRARY, GroupOption.READ_GROUP, GroupOption.MATE_CHROMOSOME,
                    GroupOption.PAIR_ORIENTATION, GroupOption.SUPPLEMENTARY, GroupOption.REFERENCE_CONCORDANCE,
                    GroupOption.MOVIE, GroupOption.ZMW, GroupOption.READ_ORDER, GroupOption.LINKED, GroupOption.PHASE
            };

            for (final GroupOption option : groupOptions) {
                JCheckBoxMenuItem mi = new JCheckBoxMenuItem(option.label);
                mi.setSelected(renderOptions.getGroupByOption() == option);
                mi.addActionListener(aEvt -> {
                    IGV.getInstance().groupAlignmentTracks(option, null, null);
                });
                groupMenu.add(mi);
                group.add(mi);
            }

            JCheckBoxMenuItem tagOption = new JCheckBoxMenuItem("tag");
            tagOption.addActionListener(aEvt -> {
                String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getGroupByTag());
                if (tag != null) {
                    if (tag.trim().length() > 0) {
                        IGV.getInstance().groupAlignmentTracks(GroupOption.TAG, tag, null);
                    } else {
                        IGV.getInstance().groupAlignmentTracks(GroupOption.NONE, null, null);
                    }
                }

            });
            tagOption.setSelected(renderOptions.getGroupByOption() == GroupOption.TAG);
            groupMenu.add(tagOption);
            group.add(tagOption);

            Range oldGroupByPos = renderOptions.getGroupByPos();
            if (oldGroupByPos != null && renderOptions.getGroupByOption() == GroupOption.BASE_AT_POS) { // already sorted by the base at a position
                JCheckBoxMenuItem oldGroupByPosOption = new JCheckBoxMenuItem("base at " + oldGroupByPos.getChr() +
                        ":" + Globals.DECIMAL_FORMAT.format(1 + oldGroupByPos.getStart()));
                groupMenu.add(oldGroupByPosOption);
                oldGroupByPosOption.setSelected(true);
            }

            if (renderOptions.getGroupByOption() != GroupOption.BASE_AT_POS || oldGroupByPos == null ||
                    !oldGroupByPos.getChr().equals(chrom) || oldGroupByPos.getStart() != chromStart) { // not already sorted by this position
                JCheckBoxMenuItem newGroupByPosOption = new JCheckBoxMenuItem("base at " + chrom +
                        ":" + Globals.DECIMAL_FORMAT.format(1 + chromStart));
                newGroupByPosOption.addActionListener(aEvt -> {
                    Range groupByPos = new Range(chrom, chromStart, chromStart + 1);
                    IGV.getInstance().groupAlignmentTracks(GroupOption.BASE_AT_POS, null, groupByPos);
                });
                groupMenu.add(newGroupByPosOption);
                group.add(newGroupByPosOption);
            }

            add(groupMenu);
        }

        /**
         * Sort menu
         */
        void addSortMenuItem() {

            JMenu sortMenu = new JMenu("Sort alignments by");
            //LinkedHashMap is supposed to preserve order of insertion for iteration
            Map<String, SortOption> mappings = new LinkedHashMap<>();

            mappings.put("start location", SortOption.START);
            mappings.put("read strand", SortOption.STRAND);
            mappings.put("first-of-pair strand", SortOption.FIRST_OF_PAIR_STRAND);
            mappings.put("base", SortOption.NUCLEOTIDE);
            mappings.put("mapping quality", SortOption.QUALITY);
            mappings.put("sample", SortOption.SAMPLE);
            mappings.put("read group", SortOption.READ_GROUP);
            mappings.put("read order", SortOption.READ_ORDER);
            mappings.put("read name", SortOption.READ_NAME);

            if (dataManager.isPairedEnd()) {
                mappings.put("insert size", SortOption.INSERT_SIZE);
                mappings.put("chromosome of mate", SortOption.MATE_CHR);
            }
            // mappings.put("supplementary flag", SortOption.SUPPLEMENTARY);

            for (Map.Entry<String, SortOption> el : mappings.entrySet()) {
                JMenuItem mi = new JMenuItem(el.getKey());
                mi.addActionListener(aEvt -> sortAlignmentTracks(el.getValue(), null));
                sortMenu.add(mi);
            }

            JMenuItem tagOption = new JMenuItem("tag");
            tagOption.addActionListener(aEvt -> {
                String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getSortByTag());
                if (tag != null && tag.trim().length() > 0) {
                    renderOptions.setSortByTag(tag);
                    sortAlignmentTracks(SortOption.TAG, tag);
                }
            });
            sortMenu.add(tagOption);

            add(sortMenu);
        }

        public void addFilterMenuItem() {
            JMenu filterMenu = new JMenu("Filter alignments by");
            JMenuItem mi = new JMenuItem("mapping quality");
            mi.addActionListener(aEvt -> {
                // TODO -- use current value for default
                String defString = PreferencesManager.getPreferences().get(SAM_QUALITY_THRESHOLD);
                if (defString == null) defString = "";
                String mqString = MessageUtils.showInputDialog("Minimum mapping quality: ", defString);
                try {
                    int mq = Integer.parseInt(mqString);
                    // TODO do something with this
                    //System.out.println(mq);
                } catch (NumberFormatException e) {
                    MessageUtils.showMessage("Mapping quality must be an integer");
                }
            });
            filterMenu.add(mi);
            add(filterMenu);
        }

        private JRadioButtonMenuItem getColorMenuItem(String label, final ColorOption option) {
            JRadioButtonMenuItem mi = new JRadioButtonMenuItem(label);
            mi.setSelected(renderOptions.getColorOption() == option);
            mi.addActionListener(aEvt -> {
                setColorOption(option);
                AlignmentTrack.this.repaint();
            });

            return mi;
        }

        void addColorByMenuItem() {
            // Change track height by attribute
            JMenu colorMenu = new JMenu("Color alignments by");

            ButtonGroup group = new ButtonGroup();

            Map<String, ColorOption> mappings = new LinkedHashMap<>();

            mappings.put("none", ColorOption.NONE);

            if (dataManager.hasYCTags()) {
                mappings.put("YC tag", ColorOption.YC_TAG);
            }

            if (dataManager.isPairedEnd()) {
                mappings.put("insert size", ColorOption.INSERT_SIZE);
                mappings.put("pair orientation", ColorOption.PAIR_ORIENTATION);
                mappings.put("insert size and pair orientation", ColorOption.UNEXPECTED_PAIR);
            }

            mappings.put("read strand", ColorOption.READ_STRAND);

            if (dataManager.isPairedEnd()) {
                mappings.put("first-of-pair strand", ColorOption.FIRST_OF_PAIR_STRAND);
            }

            mappings.put("read group", ColorOption.READ_GROUP);
            mappings.put("sample", ColorOption.SAMPLE);
            mappings.put("library", ColorOption.LIBRARY);
            mappings.put("movie", ColorOption.MOVIE);
            mappings.put("ZMW", ColorOption.ZMW);
            mappings.put("base modification", ColorOption.BASE_MODIFICATION);

            for (Map.Entry<String, ColorOption> el : mappings.entrySet()) {
                JRadioButtonMenuItem mi = getColorMenuItem(el.getKey(), el.getValue());
                colorMenu.add(mi);
                group.add(mi);
            }

            JRadioButtonMenuItem tagOption = new JRadioButtonMenuItem("tag");
            tagOption.setSelected(renderOptions.getColorOption() == ColorOption.TAG);
            tagOption.addActionListener(aEvt -> {
                setColorOption(ColorOption.TAG);
                String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getColorByTag());
                if (tag != null && tag.trim().length() > 0) {
                    setColorByTag(tag);
                    AlignmentTrack.this.repaint();
                }
            });
            colorMenu.add(tagOption);
            group.add(tagOption);

            colorMenu.add(getBisulfiteContextMenuItem(group));

            add(colorMenu);

        }

        void addPackMenuItem() {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("Re-pack alignments");
            item.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
                IGV.getInstance().packAlignmentTracks();
                AlignmentTrack.this.repaint();
            }));

            add(item);
        }

        void addCopyToClipboardItem(final TrackClickEvent te, Alignment alignment) {

            final MouseEvent me = te.getMouseEvent();
            JMenuItem item = new JMenuItem("Copy read details to clipboard");
            final ReferenceFrame frame = te.getFrame();
            if (frame == null) {
                item.setEnabled(false);
            } else {
                final double location = frame.getChromosomePosition(me.getX());

                // Change track height by attribute
                item.addActionListener(aEvt -> copyToClipboard(te, alignment, location, me.getX()));
                if (alignment == null) {
                    item.setEnabled(false);
                }
            }
            add(item);
        }


        void addViewAsPairsMenuItem() {
            final JMenuItem item = new JCheckBoxMenuItem("View as pairs");
            item.setSelected(renderOptions.isViewPairs());
            item.addActionListener(aEvt -> {
                boolean viewAsPairs = item.isSelected();
                setViewAsPairs(viewAsPairs);
            });
            item.setEnabled(dataManager.isPairedEnd());
            add(item);
        }

        void addGoToMate(final TrackClickEvent te, Alignment alignment) {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("Go to mate");
            MouseEvent e = te.getMouseEvent();

            final ReferenceFrame frame = te.getFrame();
            if (frame == null) {
                item.setEnabled(false);
            } else {
                item.addActionListener(aEvt -> gotoMate(te, alignment));
                if (alignment == null || !alignment.isPaired() || !alignment.getMate().isMapped()) {
                    item.setEnabled(false);
                }
            }
            add(item);
        }

        void showMateRegion(final TrackClickEvent te, Alignment clickedAlignment) {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("View mate region in split screen");
            MouseEvent e = te.getMouseEvent();

            final ReferenceFrame frame = te.getFrame();
            if (frame == null) {
                item.setEnabled(false);
            } else {
                double location = frame.getChromosomePosition(e.getX());

                if (clickedAlignment instanceof PairedAlignment) {
                    Alignment first = ((PairedAlignment) clickedAlignment).getFirstAlignment();
                    Alignment second = ((PairedAlignment) clickedAlignment).getSecondAlignment();
                    if (first.contains(location)) {
                        clickedAlignment = first;

                    } else if (second.contains(location)) {
                        clickedAlignment = second;

                    } else {
                        clickedAlignment = null;

                    }
                }

                final Alignment alignment = clickedAlignment;
                item.addActionListener(aEvt -> splitScreenMate(te, alignment));
                if (alignment == null || !alignment.isPaired() || !alignment.getMate().isMapped()) {
                    item.setEnabled(false);
                }
            }
            add(item);
        }

        void addClearSelectionsMenuItem() {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("Clear selections");
            item.addActionListener(aEvt -> {
                selectedReadNames.clear();
                AlignmentTrack.this.repaint();
            });
            add(item);
        }

        JMenuItem addShowAllBasesMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Show all bases");

            if (renderOptions.getColorOption() == ColorOption.BISULFITE || renderOptions.getColorOption() == ColorOption.NOMESEQ) {
                //    item.setEnabled(false);
            } else {
                item.setSelected(renderOptions.isShowAllBases());
            }
            item.addActionListener(aEvt -> {
                renderOptions.setShowAllBases(item.isSelected());
                AlignmentTrack.this.repaint();
            });
            add(item);
            return item;
        }

        void addQuickConsensusModeItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Quick consensus mode");
            item.setSelected(renderOptions.isQuickConsensusMode());

            item.addActionListener(aEvt -> {
                renderOptions.setQuickConsensusMode(item.isSelected());
                AlignmentTrack.this.repaint();
            });
            add(item);
        }

        JMenuItem addShowMismatchesMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Show mismatched bases");


            item.setSelected(renderOptions.isShowMismatches());
            item.addActionListener(aEvt -> {
                renderOptions.setShowMismatches(item.isSelected());
                AlignmentTrack.this.repaint();
            });
            add(item);
            return item;
        }


        void addInsertSizeMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Set insert size options ...");
            item.addActionListener(aEvt -> {

                InsertSizeSettingsDialog dlg = new InsertSizeSettingsDialog(IGV.getInstance().getMainFrame(), renderOptions);
                dlg.setModal(true);
                dlg.setVisible(true);
                if (!dlg.isCanceled()) {
                    renderOptions.setComputeIsizes(dlg.isComputeIsize());
                    renderOptions.setMinInsertSizePercentile(dlg.getMinPercentile());
                    renderOptions.setMaxInsertSizePercentile(dlg.getMaxPercentile());
                    if (renderOptions.computeIsizes) {
                        dataManager.updatePEStats(renderOptions);
                    }

                    renderOptions.setMinInsertSize(dlg.getMinThreshold());
                    renderOptions.setMaxInsertSize(dlg.getMaxThreshold());
                    AlignmentTrack.this.repaint();
                }
            });


            item.setEnabled(dataManager.isPairedEnd());
            add(item);
        }

        void addShadeBaseByMenuItem() {

            final JMenuItem item = new JCheckBoxMenuItem("Shade base by quality");
            item.setSelected(renderOptions.getShadeBasesOption());
            item.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
                renderOptions.setShadeBasesOption(item.isSelected());
                AlignmentTrack.this.repaint();
            }));
            add(item);
        }

        void addShowItems() {

            if (coverageTrack != null) {
                final JMenuItem item = new JCheckBoxMenuItem("Show Coverage Track");
                item.setSelected(coverageTrack.isVisible());
                item.setEnabled(!coverageTrack.isRemoved());
                item.addActionListener(aEvt -> {
                    getCoverageTrack().setVisible(item.isSelected());
                    IGV.getInstance().repaint(Arrays.asList(coverageTrack));

                });
                add(item);
            }

            if (spliceJunctionTrack != null) {
                final JMenuItem item = new JCheckBoxMenuItem("Show Splice Junction Track");
                item.setSelected(spliceJunctionTrack.isVisible());
                item.setEnabled(!spliceJunctionTrack.isRemoved());
                item.addActionListener(aEvt -> {
                    spliceJunctionTrack.setVisible(item.isSelected());
                    IGV.getInstance().repaint(Arrays.asList(spliceJunctionTrack));

                });
                add(item);
            }

            final JMenuItem alignmentItem = new JCheckBoxMenuItem("Show Alignment Track");
            alignmentItem.setSelected(true);
            alignmentItem.addActionListener(e -> {
                AlignmentTrack.this.setVisible(alignmentItem.isSelected());
                IGV.getInstance().repaint(Arrays.asList(AlignmentTrack.this));
            });
            // Disable if this is the only visible track
            if (!((coverageTrack != null && coverageTrack.isVisible()) ||
                    (spliceJunctionTrack != null && spliceJunctionTrack.isVisible()))) {
                alignmentItem.setEnabled(false);
            }

            add(alignmentItem);
        }


        void addCopySequenceItems(final TrackClickEvent te) {

            final JMenuItem item = new JMenuItem("Copy read sequence");
            add(item);
            final Alignment alignment = getSpecficAlignment(te);
            if (alignment == null) {
                item.setEnabled(false);
                return;
            }
            final String seq = alignment.getReadSequence();
            if (seq == null) {
                item.setEnabled(false);
                return;
            }
            item.addActionListener(aEvt -> StringUtils.copyTextToClipboard(seq));

            /* Add a "Copy left clipped sequence" item if there is  left clipping. */
            int minimumBlatLength = BlatClient.MINIMUM_BLAT_LENGTH;
            int[] clipping = SAMAlignment.getClipping(alignment.getCigarString());
            if (clipping[1] > 0) {
                String lcSeq = getClippedSequence(alignment.getReadSequence(), alignment.getReadStrand(), 0, clipping[1]);
                final JMenuItem lccItem = new JMenuItem("Copy left-clipped sequence");
                add(lccItem);
                lccItem.addActionListener(aEvt -> StringUtils.copyTextToClipboard(lcSeq));
            }

            /* Add a "Copy right clipped sequence" item if there is  right clipping. */
            if (clipping[3] > 0) {
                int seqLength = seq.length();
                String rcSeq = getClippedSequence(
                        alignment.getReadSequence(),
                        alignment.getReadStrand(),
                        seqLength - clipping[3],
                        seqLength);

                final JMenuItem rccItem = new JMenuItem("Copy right-clipped sequence");
                add(rccItem);
                rccItem.addActionListener(aEvt -> StringUtils.copyTextToClipboard(rcSeq));
            }
        }


        void addBlatItem(final TrackClickEvent te) {
            // Change track height by attribute
            final JMenuItem item = new JMenuItem("BLAT read sequence");
            add(item);

            final Alignment alignment = getSpecficAlignment(te);
            if (alignment == null) {
                item.setEnabled(false);
                return;
            }

            final String seq = alignment.getReadSequence();
            if (seq == null || seq.equals("*")) {
                item.setEnabled(false);
                return;

            }

            item.addActionListener(aEvt -> {
                String blatSeq = alignment.getReadStrand() == Strand.NEGATIVE ?
                        SequenceTrack.getReverseComplement(seq) : seq;
                BlatClient.doBlatQuery(blatSeq, alignment.getReadName());
            });

        }

        void addBlatClippingItems(final TrackClickEvent te) {
            final Alignment alignment = getSpecficAlignment(te);
            if (alignment == null) {
                return;
            }

            int minimumBlatLength = BlatClient.MINIMUM_BLAT_LENGTH;
            int[] clipping = SAMAlignment.getClipping(alignment.getCigarString());

            /* Add a "BLAT left clipped sequence" item if there is significant left clipping. */
            if (clipping[1] > minimumBlatLength) {
                String lcSeq = getClippedSequence(alignment.getReadSequence(), alignment.getReadStrand(), 0, clipping[1]);
                final JMenuItem lcbItem = new JMenuItem("BLAT left-clipped sequence");
                add(lcbItem);
                lcbItem.addActionListener(aEvt ->
                        BlatClient.doBlatQuery(lcSeq, alignment.getReadName() + " - left clip")
                );
            }
            /* Add a "BLAT right clipped sequence" item if there is significant right clipping. */
            if (clipping[3] > minimumBlatLength) {

                String seq = alignment.getReadSequence();
                int seqLength = seq.length();
                String rcSeq = getClippedSequence(
                        alignment.getReadSequence(),
                        alignment.getReadStrand(),
                        seqLength - clipping[3],
                        seqLength);

                final JMenuItem rcbItem = new JMenuItem("BLAT right-clipped sequence");
                add(rcbItem);
                rcbItem.addActionListener(aEvt ->
                        BlatClient.doBlatQuery(rcSeq, alignment.getReadName() + " - right clip")
                );

            }
        }

        private String getClippedSequence(String readSequence, Strand strand, int i, int i2) {
            if (readSequence == null || readSequence.equals("*")) {
                return "*";
            }
            String seq = readSequence.substring(i, i2);
            if (strand == Strand.NEGATIVE) {
                seq = SequenceTrack.getReverseComplement(seq);
            }
            return seq;
        }

        void addExtViewItem(final TrackClickEvent te) {
            // Change track height by attribute
            final JMenuItem item = new JMenuItem("ExtView");
            add(item);

            final Alignment alignment = getAlignment(te);
            if (alignment == null) {
                item.setEnabled(false);
                return;
            }

            final String seq = alignment.getReadSequence();
            if (seq == null) {
                item.setEnabled(false);
                return;

            }

            item.addActionListener(aEvt -> ExtendViewClient.postExtendView(alignment));

        }

        /**
         * Add all menu items that link alignments by tag or readname.  These are mutually exclusive.  The
         * list includes 2 items for 10X "Loupe link-read" style views, a supplementary alignment option,
         * and linking by arbitrary tag.
         */
        void addLinkedReadItems() {
            addSeparator();
            add(linkedReadViewItem("BX"));
            add(linkedReadViewItem("MI"));

            addSeparator();
            final JCheckBoxMenuItem supplementalItem = new JCheckBoxMenuItem("Link supplementary alignments");
            supplementalItem.setSelected(isLinkedReads() && "READNAME".equals(renderOptions.getLinkByTag()));
            supplementalItem.addActionListener(aEvt -> {
                boolean linkedReads = supplementalItem.isSelected();
                setLinkByTag(linkedReads, "READNAME");
            });
            add(supplementalItem);

            String linkedTagsString = PreferencesManager.getPreferences().get(SAM_LINK_BY_TAGS);
            if (linkedTagsString != null) {
                String[] t = Globals.commaPattern.split(linkedTagsString);
                for (String tag : t) {
                    if (tag.length() > 0) {
                        add(linkedReadItem(tag));
                    }
                }
            }

            final JMenuItem linkByTagItem = new JMenuItem("Link by tag...");
            linkByTagItem.addActionListener(aEvt -> {
                String tag = MessageUtils.showInputDialog("Link by tag:");
                if (tag != null) {
                    setLinkByTag(true, tag);
                    String linkedTags = PreferencesManager.getPreferences().get(SAM_LINK_BY_TAGS);
                    if (linkedTags == null) {
                        linkedTags = tag;
                    } else {
                        linkedTags += "," + tag;
                    }
                    PreferencesManager.getPreferences().put(SAM_LINK_BY_TAGS, linkedTags);
                }
            });
            add(linkByTagItem);
        }

        private JCheckBoxMenuItem linkedReadViewItem(String tag) {
            final JCheckBoxMenuItem item = new JCheckBoxMenuItem("Linked read view (" + tag + ")");
            item.setSelected(isLinkedReadView() && tag != null && tag.equals(renderOptions.getLinkByTag()));
            item.addActionListener(aEvt -> {
                boolean linkedReads = item.isSelected();
                setLinkedReadView(linkedReads, tag);
            });
            return item;
        }

        private JCheckBoxMenuItem linkedReadItem(String tag) {
            final JCheckBoxMenuItem item = new JCheckBoxMenuItem("Link by " + tag);
            item.setSelected(!isLinkedReadView() && isLinkedReads() && tag.equals(renderOptions.getLinkByTag()));
            item.addActionListener(aEvt -> {
                boolean linkedReads = item.isSelected();
                setLinkByTag(linkedReads, tag);
            });
            return item;
        }


        private void addInsertionItems(AlignmentBlock insertion) {

            final JMenuItem item = new JMenuItem("Copy insert sequence");
            add(item);
            item.addActionListener(aEvt -> StringUtils.copyTextToClipboard(insertion.getBases().getString()));

            if (insertion.getBases() != null && insertion.getBases().length >= 10) {
                final JMenuItem blatItem = new JMenuItem("BLAT insert sequence");
                add(blatItem);
                blatItem.addActionListener(aEvt -> {
                    String blatSeq = insertion.getBases().getString();
                    BlatClient.doBlatQuery(blatSeq, "BLAT insert sequence");
                });
            }
        }

        void addThirdGenItems() {

            final JMenuItem qcItem = new JCheckBoxMenuItem("Quick consensus mode");
            qcItem.setSelected(renderOptions.isQuickConsensusMode());
            qcItem.addActionListener(aEvt -> {
                renderOptions.setQuickConsensusMode(qcItem.isSelected());
                AlignmentTrack.this.repaint();
            });

            final JMenuItem thresholdItem = new JMenuItem("Small indel threshold...");
            thresholdItem.addActionListener(evt -> UIUtilities.invokeOnEventThread(() -> {
                String sith = MessageUtils.showInputDialog("Small indel threshold: ", String.valueOf(renderOptions.getSmallIndelThreshold()));
                try {
                    renderOptions.setSmallIndelThreshold(Integer.parseInt(sith));
                    AlignmentTrack.this.repaint();
                } catch (NumberFormatException e) {
                    log.error("Error setting small indel threshold - not an integer", e);
                }
            }));
            thresholdItem.setEnabled(renderOptions.isHideSmallIndels());

            final JMenuItem item = new JCheckBoxMenuItem("Hide small indels");
            item.setSelected(renderOptions.isHideSmallIndels());
            item.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
                renderOptions.setHideSmallIndels(item.isSelected());
                thresholdItem.setEnabled(item.isSelected());
                AlignmentTrack.this.repaint();
            }));

            final JMenuItem imItem = new JCheckBoxMenuItem("Show insertion markers");
            imItem.setSelected(renderOptions.isShowInsertionMarkers());
            imItem.addActionListener(aEvt -> {
                renderOptions.setShowInsertionMarkers(imItem.isSelected());
                AlignmentTrack.this.repaint();
            });

            add(imItem);
            add(qcItem);
            add(item);
            add(thresholdItem);
        }


    }


    private AlignmentBlock getInsertion(Alignment alignment, int pixelX) {
        if (alignment != null && alignment.getInsertions() != null) {
            for (AlignmentBlock block : alignment.getInsertions()) {
                if (block.containsPixel(pixelX)) {
                    return block;
                }
            }
        }
        return null;
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("experimentType")) {
            experimentType = ExperimentType.valueOf(element.getAttribute("experimentType"));
        }

        NodeList tmp = element.getElementsByTagName("RenderOptions");
        if (tmp.getLength() > 0) {
            Element renderElement = (Element) tmp.item(0);
            renderOptions = new RenderOptions(this);
            renderOptions.unmarshalXML(renderElement, version);
        }
    }


    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        if (experimentType != null) {
            element.setAttribute("experimentType", experimentType.toString());
        }

        Element sourceElement = document.createElement("RenderOptions");
        renderOptions.marshalXML(document, sourceElement);
        element.appendChild(sourceElement);

    }

    static class InsertionMenu extends IGVPopupMenu {

        final AlignmentBlock insertion;

        InsertionMenu(AlignmentBlock insertion) {

            this.insertion = insertion;

            addCopySequenceItem();

            if (insertion.getBases() != null && insertion.getBases().length > 10) {
                addBlatItem();
            }
        }


        void addCopySequenceItem() {
            // Change track height by attribute
            final JMenuItem item = new JMenuItem("Copy insert sequence");
            add(item);
            item.addActionListener(aEvt -> StringUtils.copyTextToClipboard(insertion.getBases().getString()));
        }


        void addBlatItem() {
            // Change track height by attribute
            final JMenuItem item = new JMenuItem("BLAT insert sequence");
            add(item);
            item.addActionListener(aEvt -> {
                String blatSeq = insertion.getBases().getString();
                BlatClient.doBlatQuery(blatSeq, "BLAT insert sequence");
            });
            item.setEnabled(insertion.getBases() != null && insertion.getBases().length >= 10);
        }

        @Override
        public boolean includeStandardItems() {
            return false;
        }
    }

    public static class RenderOptions implements Cloneable, Persistable {

        public static final String NAME = "RenderOptions";

        private AlignmentTrack track;
        private Boolean shadeBasesOption;
        private Boolean shadeCenters;
        private Boolean flagUnmappedPairs;
        private Boolean showAllBases;
        private Integer minInsertSize;
        private Integer maxInsertSize;
        private ColorOption colorOption;
        private GroupOption groupByOption;
        private boolean viewPairs = false;
        private String colorByTag;
        private String groupByTag;
        private String sortByTag;
        private String linkByTag;
        private Boolean linkedReads;
        private Boolean quickConsensusMode;
        private Boolean showMismatches;
        private Boolean computeIsizes;
        private Double minInsertSizePercentile;
        private Double maxInsertSizePercentile;
        private Boolean pairedArcView;
        private Boolean flagZeroQualityAlignments;
        private Range groupByPos;
        private Boolean showInsertionMarkers;
        private Boolean hideSmallIndels;
        private Integer smallIndelThreshold;

        BisulfiteContext bisulfiteContext = BisulfiteContext.CG;
        Map<String, PEStats> peStats;

        public RenderOptions() {
        }

        RenderOptions(AlignmentTrack track) {
            //updateColorScale();
            this.track = track;
            peStats = new HashMap<>();
        }

        IGVPreferences getPreferences() {
            return this.track != null ? this.track.getPreferences() : AlignmentTrack.getPreferences(ExperimentType.OTHER);
        }

        void setShowAllBases(boolean showAllBases) {
            this.showAllBases = showAllBases;
        }

        void setShowMismatches(boolean showMismatches) {
            this.showMismatches = showMismatches;
        }

        void setMinInsertSize(int minInsertSize) {
            this.minInsertSize = minInsertSize;
            //updateColorScale();
        }

        public void setViewPairs(boolean viewPairs) {
            this.viewPairs = viewPairs;
        }

        void setComputeIsizes(boolean computeIsizes) {
            this.computeIsizes = computeIsizes;
        }


        void setMaxInsertSizePercentile(double maxInsertSizePercentile) {
            this.maxInsertSizePercentile = maxInsertSizePercentile;
        }

        void setMaxInsertSize(int maxInsertSize) {
            this.maxInsertSize = maxInsertSize;
        }

        void setMinInsertSizePercentile(double minInsertSizePercentile) {
            this.minInsertSizePercentile = minInsertSizePercentile;
        }

        void setColorByTag(String colorByTag) {
            this.colorByTag = colorByTag;
        }

        void setColorOption(ColorOption colorOption) {
            this.colorOption = colorOption;
        }

        void setSortByTag(String sortByTag) {
            this.sortByTag = sortByTag;
        }

        void setGroupByTag(String groupByTag) {
            this.groupByTag = groupByTag;
        }

        void setGroupByPos(Range groupByPos) {
            this.groupByPos = groupByPos;
        }

        void setLinkByTag(String linkByTag) {
            this.linkByTag = linkByTag;
        }

        void setQuickConsensusMode(boolean quickConsensusMode) {
            this.quickConsensusMode = quickConsensusMode;
        }

        public void setGroupByOption(GroupOption groupByOption) {
            this.groupByOption = (groupByOption == null) ? GroupOption.NONE : groupByOption;
        }

        void setShadeBasesOption(boolean shadeBasesOption) {
            this.shadeBasesOption = shadeBasesOption;
        }

        void setLinkedReads(boolean linkedReads) {
            this.linkedReads = linkedReads;
        }

        public void setShowInsertionMarkers(boolean drawInsertionIntervals) {
            this.showInsertionMarkers = drawInsertionIntervals;
        }

        public void setHideSmallIndels(boolean hideSmallIndels) {
            this.hideSmallIndels = hideSmallIndels;
        }

        public void setSmallIndelThreshold(int smallIndelThreshold) {
            this.smallIndelThreshold = smallIndelThreshold;
        }

        // getters
        public int getMinInsertSize() {
            return minInsertSize == null ? getPreferences().getAsInt(SAM_MIN_INSERT_SIZE_THRESHOLD) : minInsertSize;
        }

        public int getMaxInsertSize() {
            return maxInsertSize == null ? getPreferences().getAsInt(SAM_MAX_INSERT_SIZE_THRESHOLD) : maxInsertSize;
        }

        public boolean isFlagUnmappedPairs() {
            return flagUnmappedPairs == null ? getPreferences().getAsBoolean(SAM_FLAG_UNMAPPED_PAIR) : flagUnmappedPairs;
        }

        public boolean getShadeBasesOption() {
            return shadeBasesOption == null ? getPreferences().getAsBoolean(SAM_SHADE_BASES) : shadeBasesOption;
        }

        public boolean isShowMismatches() {
            return showMismatches == null ? getPreferences().getAsBoolean(SAM_SHOW_MISMATCHES) : showMismatches;
        }

        public boolean isShowAllBases() {
            return showAllBases == null ? getPreferences().getAsBoolean(SAM_SHOW_ALL_BASES) : showAllBases;
        }

        public boolean isShadeCenters() {
            return shadeCenters == null ? getPreferences().getAsBoolean(SAM_SHADE_CENTER) : shadeCenters;
        }

        boolean isShowInsertionMarkers() {
            return showInsertionMarkers == null ? getPreferences().getAsBoolean(SAM_SHOW_INSERTION_MARKERS) : showInsertionMarkers;
        }

        public boolean isFlagZeroQualityAlignments() {
            return flagZeroQualityAlignments == null ? getPreferences().getAsBoolean(SAM_FLAG_ZERO_QUALITY) : flagZeroQualityAlignments;
        }

        public boolean isViewPairs() {
            return viewPairs;
        }

        public boolean isComputeIsizes() {
            return computeIsizes == null ? getPreferences().getAsBoolean(SAM_COMPUTE_ISIZES) : computeIsizes;
        }

        public double getMinInsertSizePercentile() {
            return minInsertSizePercentile == null ? getPreferences().getAsFloat(SAM_MIN_INSERT_SIZE_PERCENTILE) : minInsertSizePercentile;
        }

        public double getMaxInsertSizePercentile() {
            return maxInsertSizePercentile == null ? getPreferences().getAsFloat(SAM_MAX_INSERT_SIZE_PERCENTILE) : maxInsertSizePercentile;
        }

        public ColorOption getColorOption() {
            return colorOption == null ?
                    CollUtils.valueOf(ColorOption.class, getPreferences().get(SAM_COLOR_BY), ColorOption.NONE) :
                    colorOption;
        }

        public String getColorByTag() {
            return colorByTag == null ? getPreferences().get(SAM_COLOR_BY_TAG) : colorByTag;
        }

        String getSortByTag() {
            return sortByTag == null ? getPreferences().get(SAM_SORT_BY_TAG) : sortByTag;
        }

        public String getGroupByTag() {
            return groupByTag == null ? getPreferences().get(SAM_GROUP_BY_TAG) : groupByTag;
        }

        public Range getGroupByPos() {
            if (groupByPos == null) {
                String pos = getPreferences().get(SAM_GROUP_BY_POS);
                if (pos != null) {
                    String[] posParts = pos.split(" ");
                    if (posParts.length != 2) {
                        groupByPos = null;
                    } else {
                        int posChromStart = Integer.parseInt(posParts[1]);
                        groupByPos = new Range(posParts[0], posChromStart, posChromStart + 1);
                    }
                }
            }
            return groupByPos;
        }

        public String getLinkByTag() {
            return linkByTag;
        }

        public GroupOption getGroupByOption() {
            GroupOption gbo = groupByOption;
            // Interpret null as the default option.
            gbo = (gbo == null) ?
                    CollUtils.valueOf(GroupOption.class, getPreferences().get(SAM_GROUP_OPTION), GroupOption.NONE) :
                    gbo;
            // Add a second check for null in case defaultValues.groupByOption == null
            gbo = (gbo == null) ? GroupOption.NONE : gbo;

            return gbo;
        }

        public boolean isLinkedReads() {
            return linkedReads != null && linkedReads;
        }

        public boolean isQuickConsensusMode() {
            return quickConsensusMode == null ? getPreferences().getAsBoolean(SAM_QUICK_CONSENSUS_MODE) : quickConsensusMode;
        }

        public boolean isHideSmallIndels() {
            return hideSmallIndels == null ? getPreferences().getAsBoolean(SAM_HIDE_SMALL_INDEL) : hideSmallIndels;
        }

        public int getSmallIndelThreshold() {
            return smallIndelThreshold == null ? getPreferences().getAsInt(SAM_SMALL_INDEL_BP_THRESHOLD) : smallIndelThreshold;
        }


        @Override
        public void marshalXML(Document document, Element element) {

            if (shadeBasesOption != null) {
                element.setAttribute("shadeBasesOption", shadeBasesOption.toString());
            }
            if (shadeCenters != null) {
                element.setAttribute("shadeCenters", shadeCenters.toString());
            }
            if (flagUnmappedPairs != null) {
                element.setAttribute("flagUnmappedPairs", flagUnmappedPairs.toString());
            }
            if (showAllBases != null) {
                element.setAttribute("showAllBases", showAllBases.toString());
            }
            if (minInsertSize != null) {
                element.setAttribute("minInsertSize", minInsertSize.toString());
            }
            if (maxInsertSize != null) {
                element.setAttribute("maxInsertSize", maxInsertSize.toString());
            }
            if (colorOption != null) {
                element.setAttribute("colorOption", colorOption.toString());
            }
            if (groupByOption != null) {
                element.setAttribute("groupByOption", groupByOption.toString());
            }
            if (viewPairs != false) {
                element.setAttribute("viewPairs", Boolean.toString(viewPairs));
            }
            if (colorByTag != null) {
                element.setAttribute("colorByTag", colorByTag);
            }
            if (groupByTag != null) {
                element.setAttribute("groupByTag", groupByTag);
            }
            if (sortByTag != null) {
                element.setAttribute("sortByTag", sortByTag);
            }
            if (linkByTag != null) {
                element.setAttribute("linkByTag", linkByTag);
            }
            if (linkedReads != null) {
                element.setAttribute("linkedReads", linkedReads.toString());
            }
            if (quickConsensusMode != null) {
                element.setAttribute("quickConsensusMode", quickConsensusMode.toString());
            }
            if (showMismatches != null) {
                element.setAttribute("showMismatches", showMismatches.toString());
            }
            if (computeIsizes != null) {
                element.setAttribute("computeIsizes", computeIsizes.toString());
            }
            if (minInsertSizePercentile != null) {
                element.setAttribute("minInsertSizePercentile", minInsertSizePercentile.toString());
            }
            if (maxInsertSizePercentile != null) {
                element.setAttribute("maxInsertSizePercentile", maxInsertSizePercentile.toString());
            }
            if (pairedArcView != null) {
                element.setAttribute("pairedArcView", pairedArcView.toString());
            }
            if (flagZeroQualityAlignments != null) {
                element.setAttribute("flagZeroQualityAlignments", flagZeroQualityAlignments.toString());
            }
            if (groupByPos != null) {
                element.setAttribute("groupByPos", groupByPos.toString());
            }
            if (hideSmallIndels != null) {
                element.setAttribute("hideSmallIndels", hideSmallIndels.toString());
            }
            if (smallIndelThreshold != null) {
                element.setAttribute("smallIndelThreshold", smallIndelThreshold.toString());
            }
            if (showInsertionMarkers != null) {
                element.setAttribute("showInsertionMarkers", showInsertionMarkers.toString());
            }
        }


        @Override
        public void unmarshalXML(Element element, Integer version) {
            if (element.hasAttribute("shadeBasesOption")) {
                String v = element.getAttribute("shadeBasesOption");
                if (v != null) {
                    shadeBasesOption = v.equalsIgnoreCase("quality") || v.equalsIgnoreCase("true");
                }
            }
            if (element.hasAttribute("shadeCenters")) {
                shadeCenters = Boolean.parseBoolean(element.getAttribute("shadeCenters"));
            }
            if (element.hasAttribute("showAllBases")) {
                showAllBases = Boolean.parseBoolean(element.getAttribute("showAllBases"));
            }
            if (element.hasAttribute("flagUnmappedPairs")) {
                flagUnmappedPairs = Boolean.parseBoolean(element.getAttribute("flagUnmappedPairs"));
            }

            if (element.hasAttribute("minInsertSize")) {
                minInsertSize = Integer.parseInt(element.getAttribute("minInsertSize"));
            }
            if (element.hasAttribute("maxInsertSize")) {
                maxInsertSize = Integer.parseInt(element.getAttribute("maxInsertSize"));
            }
            if (element.hasAttribute("colorOption")) {
                colorOption = ColorOption.valueOf(element.getAttribute("colorOption"));
            }
            if (element.hasAttribute("groupByOption")) {
                groupByOption = GroupOption.valueOf(element.getAttribute("groupByOption"));
            }
            if (element.hasAttribute("viewPairs")) {
                viewPairs = Boolean.parseBoolean(element.getAttribute("viewPairs"));
            }
            if (element.hasAttribute("colorByTag")) {
                colorByTag = element.getAttribute("colorByTag");
            }
            if (element.hasAttribute("groupByTag")) {
                groupByTag = element.getAttribute("groupByTag");
            }
            if (element.hasAttribute("sortByTag")) {
                sortByTag = element.getAttribute("sortByTag");
            }
            if (element.hasAttribute("linkByTag")) {
                linkByTag = element.getAttribute("linkByTag");
            }
            if (element.hasAttribute("linkedReads")) {
                linkedReads = Boolean.parseBoolean(element.getAttribute("linkedReads"));
            }
            if (element.hasAttribute("quickConsensusMode")) {
                quickConsensusMode = Boolean.parseBoolean(element.getAttribute("quickConsensusMode"));
            }
            if (element.hasAttribute("showMismatches")) {
                showMismatches = Boolean.parseBoolean(element.getAttribute("showMismatches"));
            }
            if (element.hasAttribute("computeIsizes")) {
                computeIsizes = Boolean.parseBoolean(element.getAttribute("computeIsizes"));
            }
            if (element.hasAttribute("minInsertSizePercentile")) {
                minInsertSizePercentile = Double.parseDouble(element.getAttribute("minInsertSizePercentile"));
            }
            if (element.hasAttribute("maxInsertSizePercentile")) {
                maxInsertSizePercentile = Double.parseDouble(element.getAttribute("maxInsertSizePercentile"));
            }
            if (element.hasAttribute("pairedArcView")) {
                pairedArcView = Boolean.parseBoolean(element.getAttribute("pairedArcView"));
            }
            if (element.hasAttribute("flagZeroQualityAlignments")) {
                flagZeroQualityAlignments = Boolean.parseBoolean(element.getAttribute("flagZeroQualityAlignments"));
            }
            if (element.hasAttribute("groupByPos")) {
                groupByPos = Range.fromString(element.getAttribute("groupByPos"));
            }
            if (element.hasAttribute("hideSmallIndels")) {
                hideSmallIndels = Boolean.parseBoolean(element.getAttribute("hideSmallIndels"));
            }
            if (element.hasAttribute("smallIndelThreshold")) {
                smallIndelThreshold = Integer.parseInt(element.getAttribute("smallIndelThreshold"));
            }
            if (element.hasAttribute("showInsertionMarkers")) {
                showInsertionMarkers = Boolean.parseBoolean(element.getAttribute("showInsertionMarkers"));
            }
        }
    }


}
