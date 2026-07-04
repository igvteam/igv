package org.igv.alignment;


import org.igv.Globals;
import org.igv.event.*;
import org.igv.feature.FeatureUtils;
import org.igv.feature.Range;
import org.igv.feature.genome.Genome;
import org.igv.circview.CircularViewUtilities;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.GraphicUtils;
import org.igv.track.*;
import org.igv.ui.FontManager;
import org.igv.ui.IGV;
import org.igv.ui.color.ColorTable;
import org.igv.ui.color.ColorUtilities;
import org.igv.ui.color.PaletteColorTable;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.MessageUtils;
import org.igv.util.ResourceLocator;
import org.igv.util.StringUtils;
import org.igv.util.blat.BlatClient;
import org.json.JSONObject;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.util.*;
import java.util.List;

import static org.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */

public class AlignmentTrack extends AbstractTrack implements IGVEventObserver {

    private static final Logger log = LogManager.getLogger(AlignmentTrack.class);

    static final int LEGACY_SQUISHED_HEIGHT = 2;

    // Alignment colors
    static final Color DEFAULT_ALIGNMENT_COLOR = new Color(185, 185, 185); //200, 200, 200);

    public enum ColorOption {
        INSERT_SIZE,
        READ_STRAND,
        FIRST_OF_PAIR_STRAND,
        PAIR_ORIENTATION,
        READ_ORDER,
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
        SPLIT,
        BASE_MODIFICATION,
        BASE_MODIFICATION_2COLOR,
        SMRT_SUBREAD_IPD,
        SMRT_SUBREAD_PW,
        SMRT_CCS_FWD_IPD,
        SMRT_CCS_FWD_PW,
        SMRT_CCS_REV_IPD,
        SMRT_CCS_REV_PW;

        public boolean isBaseMod() {
            return this == BASE_MODIFICATION || this == BASE_MODIFICATION_2COLOR;
        }

        public boolean isSMRTKinetics() {
            switch (this) {
                case SMRT_SUBREAD_IPD:
                case SMRT_SUBREAD_PW:
                case SMRT_CCS_FWD_IPD:
                case SMRT_CCS_REV_IPD:
                case SMRT_CCS_FWD_PW:
                case SMRT_CCS_REV_PW:
                    return true;
                default:
                    return false;
            }
        }
    }

    public enum DuplicatesOption {
        FILTER("filter duplicates", true),
        SHOW("show duplicates", false),
        TEXTURE("texture duplicates", false);

        public final String label;
        public final boolean filtered;

        DuplicatesOption(String label, Boolean filtered) {
            this.label = label;
            this.filtered = filtered;
        }
    }

    public enum ShadeAlignmentsOption {
        NONE("none"),
        MAPPING_QUALITY_HIGH("mapping quality high"),
        MAPPING_QUALITY_LOW("mapping quality low");

        public final String label;

        ShadeAlignmentsOption(String label) {
            this.label = label;
        }
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
        CHIMERIC("chimeric"),
        SUPPLEMENTARY("supplementary flag"),
        BASE_AT_POS("base at position"),
        INSERTION_AT_POS("insertion at position", true),
        MOVIE("movie"),
        ZMW("ZMW"),
        CLUSTER("cluster"),
        READ_ORDER("read order"),
        LINKED("linked"),
        PHASE("phase"),
        REFERENCE_CONCORDANCE("reference concordance"),
        MAPPING_QUALITY("mapping quality"),
        SELECTED("selected"),
        DUPLICATE("duplicate flag");

        public final String label;
        public final boolean reverse;

        GroupOption(String label) {
            this.label = label;
            this.reverse = false;
        }

        GroupOption(String label, boolean reverse) {
            this.label = label;
            this.reverse = reverse;
        }

    }

    enum OrientationType {
        RR, LL, RL, LR, UNKNOWN
    }

    private static final int GROUP_LABEL_HEIGHT = 10;
    private static final int GROUP_MARGIN = 5;
    private static final int TOP_MARGIN = 20;
    private static final int DS_MARGIN_0 = 2;
    private static final int DOWNSAMPLED_ROW_HEIGHT = 5;

    public enum BisulfiteContext {
        CG("CG", new byte[]{}, new byte[]{'G'}),
        CHH("CHH", new byte[]{}, new byte[]{'H', 'H'}),
        CHG("CHG", new byte[]{}, new byte[]{'H', 'G'}),
        HCG("HCG", new byte[]{'H'}, new byte[]{'G'}),
        GCH("GCH", new byte[]{'G'}, new byte[]{'H'}),
        WCG("WCG", new byte[]{'W'}, new byte[]{'G'}),
        NONE("None", null, null);

        private final String label;
        private final byte[] preContext;
        private final byte[] postContext;

        /**
         * @param contextb      The residue in the context string (IUPAC)
         * @param referenceBase The reference sequence (already checked that offsetidx is within bounds)
         * @param readBase      The read sequence (already checked that offsetidx is within bounds)
         */
        private static boolean positionMatchesContext(byte contextb, final byte referenceBase, final byte readBase) {
            boolean matchesContext = AlignmentUtils.compareBases(contextb, referenceBase);
            if (!matchesContext) {
                return false; // Don't need to check any further
            }

            // For the read, we have to handle C separately
            boolean matchesReadContext = AlignmentUtils.compareBases(contextb, readBase);
            if (AlignmentUtils.compareBases((byte) 'T', readBase)) {
                matchesReadContext |= AlignmentUtils.compareBases(contextb, (byte) 'C');
            }

            return matchesReadContext;
        }

        public BisulfiteContext getMatchingBisulfiteContext(final byte[] reference, final ByteSubarray read, final int idx) {
            boolean matchesContext = true;

            // First do the "post" context
            int minLen = Math.min(reference.length, read.length);
            if ((idx + postContext.length) >= minLen) {
                matchesContext = false;
            } else {
                // Cut short whenever we don't match
                for (int posti = 0; matchesContext && (posti < postContext.length); posti++) {
                    byte contextb = postContext[posti];
                    int offsetidx = idx + 1 + posti;

                    matchesContext &= positionMatchesContext(contextb, reference[offsetidx], read.getByte(offsetidx));
                }
            }

            // Now do the pre context
            if ((idx - preContext.length) < 0) {
                matchesContext = false;
            } else {
                // Cut short whenever we don't match
                for (int prei = 0; matchesContext && (prei < preContext.length); prei++) {
                    byte contextb = preContext[prei];
                    int offsetidx = idx - (preContext.length - prei);

                    matchesContext &= positionMatchesContext(contextb, reference[offsetidx], read.getByte(offsetidx));
                }
            }

            return (matchesContext) ? this : null;
        }

        public String getLabel() {
            return label;
        }

        BisulfiteContext(String label, byte[] preContext, byte[] postContext) {
            this.label = label;
            this.preContext = preContext;
            this.postContext = postContext;
        }
    }

    public static boolean isBisulfiteColorType(ColorOption o) {
        return (o.equals(ColorOption.BISULFITE) || o.equals(ColorOption.NOMESEQ));
    }

    private final AlignmentDataManager dataManager;
    private final SequenceTrack sequenceTrack;
    private final CoverageTrack coverageTrack;
    private final SpliceJunctionTrack spliceJunctionTrack;

    private final Genome genome;
    ExperimentType experimentType;
    private final AlignmentRenderer renderer;
    RenderOptions renderOptions;

    private boolean removed = false;
    private boolean showGroupLine;
    private int collapsedHeight = 9;
    private int squishedHeight = 2;
    private final int minHeight = 50;

    private Rectangle downsampleRect;
    private ColorTable readNamePalette;
    private final HashMap<String, Color> selectedReadNames = new HashMap<>();

    /**
     * Create a new alignment track
     *
     * @param locator
     * @param dataManager
     * @param genome
     */
    public AlignmentTrack(ResourceLocator locator, AlignmentDataManager dataManager, Genome genome) {
        super(locator);

        final String baseName = locator.getTrackName();
        this.setName(baseName);
        this.rowHeight = 14;
        this.dataManager = dataManager;
        this.genome = genome;
        this.renderer = new AlignmentRenderer(this);
        this.renderOptions = new RenderOptions(this);
        setColor(DEFAULT_ALIGNMENT_COLOR);
        dataManager.setAlignmentTrack(this);
        dataManager.subscribe(this);

        IGVPreferences prefs = getPreferences();
        setHeight(PreferencesManager.getPreferences().getAsInt(ALIGNMENT_TRACK_HEIGHT));
        minimumHeight = 50;
        showGroupLine = prefs.getAsBoolean(SAM_SHOW_GROUP_SEPARATOR);
        try {
            setDisplayMode(DisplayMode.valueOf(prefs.get(SAM_DISPLAY_MODE).toUpperCase()));
        } catch (Exception e) {
            setDisplayMode(DisplayMode.EXPANDED);
        }

        // Optional sequence track (not common)
        if (prefs.getAsBoolean(SAM_SHOW_REF_SEQ)) {
            sequenceTrack = new SequenceTrack("Reference sequence");
            sequenceTrack.setHeight(14);
        } else {
            sequenceTrack = null;
        }

        // Coverage track
        this.coverageTrack = new CoverageTrack(locator, baseName + " Coverage", this, genome);
        this.coverageTrack.setDataManager(dataManager);
        dataManager.setCoverageTrack(this.coverageTrack);

        // Splice junction track
        SpliceJunctionTrack spliceJunctionTrack = new SpliceJunctionTrack(locator,
                baseName + " Junctions", dataManager, this, SpliceJunctionTrack.StrandOption.BOTH);
        spliceJunctionTrack.setHeight(60);
        this.spliceJunctionTrack = spliceJunctionTrack;

        if (renderOptions.getColorOption() == ColorOption.BISULFITE) {
            setExperimentType(ExperimentType.BISULFITE);
        }
        readNamePalette = new PaletteColorTable(ColorUtilities.getDefaultPalette());

        dataManager.setViewAsPairs(prefs.getAsBoolean(SAM_DISPLAY_PAIRED), renderOptions, getDisplayMode());

        IGVEventBus.getInstance().subscribe(FrameManager.ChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(AlignmentTrackEvent.class, this);
        IGVEventBus.getInstance().subscribe(DataLoadedEvent.class, this);
        IGVEventBus.getInstance().subscribe(ViewChange.class, this);
    }

    public void init() {
        if (experimentType == null || experimentType == ExperimentType.UNKOWN) {
            ExperimentType type = dataManager.inferType();
            setExperimentType(type);
        }
    }

    @Override
    public TrackType getType() {
        return TrackType.alignment;
    }

    @Override
    public void receiveEvent(IGVEvent event) {

        if (event instanceof FrameManager.ChangeEvent) {
            // Trim insertionInterval map to current frames

        } else if (event instanceof AlignmentTrackEvent e) {
            switch (e.type()) {
                case ALLELE_THRESHOLD -> dataManager.alleleThresholdChanged();
                case RELOAD -> {
                    clearCaches();
                    repaint();
                }
                case REFRESH -> repaint();
            }
        } else if (event instanceof DataLoadedEvent dataLoaded) {
            if (getPreferences().getAsBoolean(SAM_AUTO_SORT)) {
                sortRows(dataLoaded.referenceFrame());
            }
        } else if (event instanceof ViewChange viewChange) {
            if (viewChange.type == ViewChange.Type.LocusChange && !viewChange.panning) {
                if (getDisplayMode() == DisplayMode.FULL) {
                    packAlignments();
                }
                // Don't autosort on completion of a track pan (drag)
                if (getPreferences().getAsBoolean(SAM_AUTO_SORT) && !viewChange.panning) {
                    sortRows(viewChange.referenceFrame);
                }
            }
        }
    }

    void setExperimentType(ExperimentType type) {

        if (type != experimentType) {

            experimentType = type;
            boolean revalidate = false;

            boolean showJunction = getPreferences(type).getAsBoolean(Constants.SAM_SHOW_JUNCTION_TRACK);
            if (showJunction != spliceJunctionTrack.isVisible()) {
                spliceJunctionTrack.setVisible(showJunction);
                revalidate = true;
            }

            boolean showCoverage = getPreferences(type).getAsBoolean(SAM_SHOW_COV_TRACK);
            if (showCoverage != coverageTrack.isVisible()) {
                coverageTrack.setVisible(showCoverage);
                revalidate = true;
            }

            boolean showAlignments = getPreferences(type).getAsBoolean(SAM_SHOW_ALIGNMENT_TRACK);
            if (showAlignments != isVisible()) {
                setVisible(showAlignments);
                revalidate = true;
            }

            if (IGV.hasInstance()) {
                if (revalidate) {
                    IGV.getInstance().revalidateTrackPanels();
                } else {
                    this.repaint();
                }
            }

            //ExperimentTypeChangeEvent event = new ExperimentTypeChangeEvent(this, experimentType);
            //IGVEventBus.getInstance().post(event);
        }
    }

    public ExperimentType getExperimentType() {
        return experimentType;
    }

    public AlignmentDataManager getDataManager() {
        return dataManager;
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public SpliceJunctionTrack getSpliceJunctionTrack() {
        return spliceJunctionTrack;
    }

    public RenderOptions getRenderOptions() {
        return renderOptions;
    }

    public HashMap<String, Color> getSelectedReadNames() {
        return selectedReadNames;
    }

    public ColorTable getReadNamePalette() {
        return readNamePalette;
    }

    @Override
    public List<Component> getPopupMenuItems(TrackClickEvent te) {
        return AlignmentTrackMenuHelper.getMenuItems(this, te);
    }

    @Override
    public void setDisplayMode(DisplayMode mode) {
        // Transitioning to or from FULL requires repacking
        boolean repack = (getDisplayMode() == DisplayMode.FULL || mode == DisplayMode.FULL);

        // Legacy "squished" mode -- an expanded mode with reduced row height
        if(mode == DisplayMode.SQUISHED) {
            setRowHeight(LEGACY_SQUISHED_HEIGHT);
            mode = DisplayMode.EXPANDED;
        }

        super.setDisplayMode(mode);
        if (repack) {
            packAlignments();
        }
    }


    @Override
    public int getContentHeight() {

        int nGroups = dataManager.getMaxGroupCount();
        int intRowHeight = Math.max(1, rowHeight);  // This is what is used to draw
        int h = Math.max(minHeight, getNLevels() * intRowHeight + nGroups * GROUP_MARGIN + TOP_MARGIN
                + DS_MARGIN_0 + DOWNSAMPLED_ROW_HEIGHT);
        return Math.max(minimumHeight, h);
    }

    @Override
    public int getRowHeight() {
        final DisplayMode displayMode = getDisplayMode();
        if (displayMode == DisplayMode.EXPANDED || displayMode == DisplayMode.FULL) {
            return rowHeight;
        } else if (displayMode == DisplayMode.COLLAPSED) {
            return collapsedHeight;
        } else {
            return squishedHeight;
        }
    }

    @Override
    public void minimizeHeight() {
        setRowHeight(1);
        int newHeight = Math.max(getContentHeight(), getMinimumHeight());
        setHeight(Math.min(newHeight, getHeight()));
    }

    @Override
    public int getNumRows() {
        return getNLevels();
    }

    private int getNLevels() {
        return dataManager.getNLevels();
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {

        double extent = frame.getEnd() - frame.getOrigin();
        if (frame.getChrName().equals(Globals.CHR_ALL) || extent > getVisibilityWindow()) {
            return true;   // Nothing to paint
        } else {
            return dataManager.isLoaded(frame);
        }
    }


    @Override
    public void load(ReferenceFrame referenceFrame) {
        if (log.isDebugEnabled()) {
            log.debug("Reading - thread: " + Thread.currentThread().getName());
        }
        dataManager.load(referenceFrame, renderOptions, getDisplayMode(), true);
        sortRows(referenceFrame);
    }

    @Override
    public int getVisibilityWindow() {
        return (int) dataManager.getVisibilityWindow();
    }


    /**
     * Draw that portion of the alignments track overlapping visibleRect.  This is the subset of the track visible
     * in the current viewport.
     *
     * @param context
     */


    public void render(RenderContext context) {

        int viewWindowSize = context.getReferenceFrame().getCurrentRange().getLength();
        if (viewWindowSize > getVisibilityWindow()) {
            return;
        }

        Rectangle trackRect = context.getTrackRectangle();
        Rectangle clipBounds = context.getClipBounds();
        context.getGraphics2D("LABEL").setFont(FontManager.getFont(GROUP_LABEL_HEIGHT));

        int seqHeight = sequenceTrack == null ? 0 : sequenceTrack.getHeight();
        int yOffset = seqHeight;
        if (clipBounds.y < seqHeight) {
            if (seqHeight > 0) {
                Rectangle seqRect = new Rectangle(0, 0, trackRect.width, seqHeight);
                sequenceTrack.render(context);
            }
        }

        // Top gap.
        boolean downsampled = false;
        if (clipBounds.y < yOffset + DS_MARGIN_0 + DOWNSAMPLED_ROW_HEIGHT) {
            Rectangle dsRect = new Rectangle(trackRect);
            dsRect.y = yOffset + DS_MARGIN_0;
            dsRect.height = DOWNSAMPLED_ROW_HEIGHT;
            downsampled = renderDownsampledIntervals(context, dsRect);
        }
        if (downsampled) {
            yOffset += DS_MARGIN_0 + DOWNSAMPLED_ROW_HEIGHT;
            this.downsampleRect = downsampleRect;
        } else {
            this.downsampleRect = null;
        }
        yOffset += DS_MARGIN_0;

        Rectangle alignmentsRect = new Rectangle(trackRect);
        alignmentsRect.y = yOffset;
        alignmentsRect.height -= yOffset;
        renderAlignments(context, alignmentsRect);
    }

    private boolean renderDownsampledIntervals(RenderContext context, Rectangle downsampleRect) {

        // Might be offscreen

        final AlignmentInterval loadedInterval = dataManager.getLoadedInterval(context.getReferenceFrame());
        if (loadedInterval == null) return false;

        Graphics2D g = context.getGraphic2DForColor(darkMode ? Color.white : Color.black);

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
        return intervals.size() > 0;
    }


    private void renderAlignments(RenderContext context, Rectangle alignmentsRect) {

        Rectangle clipBounds = context.getClipBounds();

        final AlignmentInterval loadedInterval = dataManager.getLoadedInterval(context.getReferenceFrame(), true);
        if (loadedInterval == null) {
            return;
        }

        final AlignmentCounts alignmentCounts = loadedInterval.getCounts();

        PackedAlignments groups = dataManager.getGroups(loadedInterval, renderOptions);
        if (groups == null) {
            //Assume we are still loading.
            return;
        }

        // Check for YC tag
        if (renderOptions.getColorOption() == null && dataManager.hasYCTags()) {
            renderOptions.setColorOption(ColorOption.YC_TAG);
        }

        Map<String, PEStats> peStats = dataManager.getPEStats();
        if (peStats != null) {
            renderOptions.peStats = peStats;
        }

        int intH = Math.max(1, rowHeight);

        // Loop through groups
        double y = alignmentsRect.y;
        Graphics2D groupBorderGraphics = context.getGraphic2DForColor(AlignmentRenderer.GROUP_DIVIDER_COLOR);
        int nGroups = groups.size();
        int groupNumber = 0;
        GroupOption groupOption = renderOptions.getGroupByOption();
        for (Map.Entry<String, List<Row>> entry : groups.entrySet()) {

            groupNumber++;
            double yGroup = y;  // Remember this for label

            // Loop through the alignment rows for this group
            // Create a snapshot to avoid ConcurrentModificationException from background loading threads
            List<Row> rows = new ArrayList<>(entry.getValue());
            for (Row row : rows) {
                if (y > clipBounds.getMaxY()) {
                    break;
                }

                if (y + intH > clipBounds.y) {
                    Rectangle rowRectangle = new Rectangle(alignmentsRect.x, (int) y, alignmentsRect.width, intH);
                    renderer.renderAlignments(row.alignments, alignmentCounts, context, rowRectangle, renderOptions);
                    row.y = y;
                    row.h = intH;
                }
                y += intH;
            }

            if (groupOption != GroupOption.NONE) {
                // Draw a subtle divider line between groups
                if (showGroupLine) {
                    if (groupNumber < nGroups) {
                        int borderY = (int) y + GROUP_MARGIN / 2;
                        GraphicUtils.drawDottedDashLine(groupBorderGraphics, alignmentsRect.x, borderY, alignmentsRect.width, borderY);
                    }
                }

                // Label the group, if there is room
                double groupHeight = rows.size() * intH;
                if (groupHeight > GROUP_LABEL_HEIGHT + 2 && !context.multiframe) {
                    String groupName = entry.getKey();
                    if (groupName.equals("SELECTED")) {
                        // Abbreviate the "SELECTED" group label to "S*" for display to conserve horizontal space in the UI.
                        groupName = "S*";
                    }
                    Graphics2D g = context.getGraphics2D("LABEL");
                    FontMetrics fm = g.getFontMetrics();
                    Rectangle2D stringBouds = fm.getStringBounds(groupName, g);
                    Rectangle rect = new Rectangle(alignmentsRect.x, (int) yGroup, (int) stringBouds.getWidth() + 10, (int) stringBouds.getHeight());
                    GraphicUtils.drawVerticallyCenteredText(groupName, 5, rect, g, false, true);
                }
            }
            y += GROUP_MARGIN;
        }

        final int bottom = alignmentsRect.y + alignmentsRect.height;
        groupBorderGraphics.drawLine(alignmentsRect.x, bottom, alignmentsRect.width, bottom);
    }

    /**
     * Render insertions at position of marker in expanded form, showing sequence.  Much of this method is just
     * a repeat of the loop for alignment rendering to compute the Y position of affected alignments.
     *
     * @param insertionMarker
     * @param context
     * @param inputRect
     */
    public void renderExpandedInsertion(InsertionMarker insertionMarker, RenderContext context, Rectangle inputRect) {

        boolean leaveMargin = rowHeight > 2;
        inputRect.y += DS_MARGIN_0 + DOWNSAMPLED_ROW_HEIGHT + DS_MARGIN_0;

        final AlignmentInterval loadedInterval = dataManager.getLoadedInterval(context.getReferenceFrame(), true);
        PackedAlignments groups = dataManager.getGroups(loadedInterval, renderOptions);
        if (groups == null) {
            //Assume we are still loading.
            return;
        }

        Rectangle clipBounds = context.getClipBounds();

        // Divide rectangle into equal height levels
        double y = inputRect.getY() - 3;
        int intH;
        if (getDisplayMode() == DisplayMode.EXPANDED) {
            intH = Math.max(1, rowHeight);
        } else if (getDisplayMode() == DisplayMode.COLLAPSED) {
            intH = collapsedHeight;
        } else {
            intH = squishedHeight;
        }

        for (Map.Entry<String, List<Row>> entry : groups.entrySet()) {
            // Loop through the alignment rows for this group
            List<Row> rows = entry.getValue();
            for (Row row : rows) {
                if ((clipBounds != null && y > clipBounds.getMaxY())) {
                    return;
                }

                assert clipBounds != null;
                if (y + intH > clipBounds.getY()) {
                    Rectangle rowRectangle = new Rectangle(inputRect.x, (int) y, inputRect.width, intH);
                    if (row.alignments != null)  // TODO -- not sure this is needed
                        BaseRenderer.drawExpandedInsertions(insertionMarker, row.alignments, context, rowRectangle, leaveMargin, renderOptions);
                    row.y = y;
                    row.h = intH;
                }
                y += intH;
            }

            y += GROUP_MARGIN;

        }
    }

    /**
     * Sort alignment rows for all reference frames.  This is called on user initiated sort events.
     *
     * @param option
     * @param tag
     * @param invertSort
     */
    public void sortRows(final SortOption option, final String tag, final boolean invertSort) {
        final List<ReferenceFrame> frames = FrameManager.getFrames();
        for (ReferenceFrame frame : frames) {
            final AlignmentInterval interval = getDataManager().getLoadedInterval(frame);
            if (interval != null) {
                final double location = frame.getCenter();
                interval.sortRows(option, location, tag, invertSort);
            }
        }
    }

    /**
     * Sort alignment rows for a specific reference frame using track render settings.  This is called on locus
     * change and data load events.
     *
     * @param referenceFrame
     */
    public void sortRows(ReferenceFrame referenceFrame) {
        SortOption option = this.renderOptions.getSortOption();
        String tag = this.renderOptions.getSortByTag();
        boolean invertSort = this.renderOptions.isInvertSorting();
        if (option != SortOption.NONE) {
            final AlignmentInterval interval = getDataManager().getLoadedInterval(referenceFrame);
            if (interval != null) {
                double location = referenceFrame.getCenter();
                interval.sortRows(option, location, tag, invertSort);
            }
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
        if ((option == GroupOption.BASE_AT_POS || option == GroupOption.INSERTION_AT_POS) && pos != null) {
            renderOptions.setGroupByPos(pos);
        }
        renderOptions.setGroupByOption(option);
        dataManager.packAlignments(renderOptions, getDisplayMode());
        repaint();
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

    public void setShadeAlignmentsOptions(ShadeAlignmentsOption option) {
        renderOptions.setShadeAlignmentsOption(option);
    }

    public void packAlignments() {
        dataManager.packAlignments(renderOptions, getDisplayMode());
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
            Alignment feature = getAlignmentAt(position, mouseY, frame);
            if (feature != null) {
                return feature.getAlignmentValueString(position, mouseX, renderOptions);
            }
        }
        return null;
    }


    Alignment getAlignmentAt(final TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();
        final ReferenceFrame frame = te.getFrame();
        return frame == null ? null : getAlignmentAt(frame.getChromosomePosition(e), e.getY(), frame);
    }

    Alignment getAlignmentAt(double position, int y, ReferenceFrame frame) {

        if (dataManager == null) {
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


    @Override
    public boolean handleDataClick(TrackClickEvent te) {

        MouseEvent e = te.getMouseEvent();

        if (Globals.IS_MAC && e.isMetaDown() || (!Globals.IS_MAC && e.isControlDown())) {
            // Selection
            Alignment alignment = this.getAlignmentAt(te);
            if (alignment != null) {
                if (selectedReadNames.containsKey(alignment.getReadName())) {
                    selectedReadNames.remove(alignment.getReadName());
                } else {
                    setSelectedAlignment(alignment);
                }
                IGV.getInstance().repaint(this); //todo check if doing this conditionally here is ok
            }
            return true;
        }

        if (IGV.getInstance().isShowDetailsOnClick()) {
            openTooltipWindow(te);
            return true;
        }

        return false;
    }

    public void setSelectedAlignment(Alignment alignment) {
        Color c = readNamePalette.get(alignment.getReadName());
        selectedReadNames.put(alignment.getReadName(), c);
    }

    private void clearCaches() {
        if (dataManager != null) dataManager.clear();
        if (spliceJunctionTrack != null) spliceJunctionTrack.clear();
    }


    public void setViewAsPairs(boolean vAP) {
        // TODO -- generalize this test to all incompatible pairings
        if (vAP && renderOptions.getGroupByOption() == GroupOption.STRAND) {
            boolean ungroup = MessageUtils.confirm("\"View as pairs\" is incompatible with \"Group by strand\". Ungroup?");
            if (ungroup) {
                renderOptions.setGroupByOption(null);
            } else {
                return;
            }
        }

        dataManager.setViewAsPairs(vAP, renderOptions, getDisplayMode());
        repaint();
    }


    public enum ExperimentType {OTHER, RNA, BISULFITE, THIRD_GEN, SBX, UNKOWN}


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

    boolean isLinkedReads() {
        return renderOptions != null && renderOptions.isLinkedReads();
    }

    /**
     * Set the view to a 10X style "linked-read view".  This option tries to achieve a view similar to the
     * 10X Loupe view described here: https://support.10xgenomics.com/genome-exome/software/visualization/latest/linked-reads
     *
     * @param linkedReads
     * @param tag
     */
    void setLinkedReadView(boolean linkedReads, String tag) {
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
            setRowHeight(2);
            setDisplayMode(DisplayMode.EXPANDED);
        }
        dataManager.packAlignments(renderOptions, getDisplayMode());
        repaint();
    }

    /**
     * Detect if we are in linked-read view (10X Loupe style view)
     */
    boolean isLinkedReadView() {
        return renderOptions != null &&
                renderOptions.isLinkedReads() &&
                renderOptions.getLinkByTag() != null &&
                renderOptions.getColorOption() == ColorOption.TAG &&
                renderOptions.getColorByTag() != null;
    }

    void undoLinkedReadView() {
        renderOptions.setLinkByTag(null);
        renderOptions.setColorOption(ColorOption.NONE);
        renderOptions.setColorByTag(null);
        renderOptions.setGroupByOption(GroupOption.NONE);
        renderOptions.setGroupByTag(null);
        showGroupLine = true;
        setDisplayMode(DisplayMode.EXPANDED);
    }

    void sendPairsToCircularView(TrackClickEvent e) {

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
            CircularViewUtilities.addAlignments(inView, AlignmentTrack.this.getName(), chordColor);
        }
    }

    void sendSplitToCircularView(TrackClickEvent e) {

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
            CircularViewUtilities.addAlignments(inView, AlignmentTrack.this.getName(), chordColor);
        }
    }


    AlignmentBlock getInsertion(Alignment alignment, int pixelX) {
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
    public void unmarshalJSON(JSONObject jsonObject) {
        super.unmarshalJSON(jsonObject);

        if (jsonObject.has("experimentType")) {
            experimentType = ExperimentType.valueOf(jsonObject.getString("experimentType"));
        }

        renderOptions = new RenderOptions(this);
        renderOptions.unmarshalJSON(jsonObject);
    }

    @Override
    public void marshalJSON(JSONObject jsonObject) {
        super.marshalJSON(jsonObject);

        if (experimentType != null) {
            jsonObject.put("experimentType", experimentType.toString());
        }

        renderOptions.marshalJSON(jsonObject);
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


}
