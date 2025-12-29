package org.igv.sam;


import org.igv.Globals;
import org.igv.event.*;
import org.igv.feature.FeatureUtils;
import org.igv.feature.Range;
import org.igv.feature.genome.Genome;
import org.igv.jbrowse.CircularViewUtilities;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.renderer.GraphicUtils;
import org.igv.sam.mods.BaseModficationFilter;
import org.igv.session.Persistable;
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
import org.igv.util.collections.CollUtils;
import org.w3c.dom.Document;
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
    private static final int DOWNSAMPLED_ROW_HEIGHT = 10;
    private static final int INSERTION_ROW_HEIGHT = 9;

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
    private ExperimentType experimentType;
    private final AlignmentRenderer renderer;
    RenderOptions renderOptions;

    private boolean removed = false;
    private boolean showGroupLine;
    private int expandedHeight = 14;
    private int collapsedHeight = 9;
    private final int maxSquishedHeight = 5;
    private int squishedHeight = maxSquishedHeight;
    private final int minHeight = 50;

    private Rectangle alignmentsRect;
    private Rectangle downsampleRect;
    private Rectangle insertionRect;
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
        this.dataManager = dataManager;
        this.genome = genome;
        this.renderer = new AlignmentRenderer(this);
        this.renderOptions = new RenderOptions(this);
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
    public boolean isAlignment() {
        return true;
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
            sortRows(dataLoaded.referenceFrame());
        } else if (event instanceof ViewChange viewChange) {
            if (viewChange.type == ViewChange.Type.LocusChange && !viewChange.panning) {
                if (getDisplayMode() == DisplayMode.FULL) {
                    packAlignments();
                }
                // Don't autosort on completion of a track pan (drag)
                if (!viewChange.fromPanning) {
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
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {
        return new AlignmentTrackMenu(this, te);
    }

    @Override
    public void setDisplayMode(DisplayMode mode) {
        boolean repack = (getDisplayMode() == DisplayMode.FULL || mode == DisplayMode.FULL);
        super.setDisplayMode(mode);
        if (repack) {
            packAlignments();
        }
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
                + DS_MARGIN_0 + DOWNSAMPLED_ROW_HEIGHT);
        return Math.max(minimumHeight, h);
    }

    private int getRowHeight() {
        final DisplayMode displayMode = getDisplayMode();
        if (displayMode == DisplayMode.EXPANDED || displayMode == DisplayMode.FULL) {
            return expandedHeight;
        } else if (displayMode == DisplayMode.COLLAPSED) {
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

    public void render(RenderContext context, Rectangle rect) {

        int viewWindowSize = context.getReferenceFrame().getCurrentRange().getLength();
        if (viewWindowSize > getVisibilityWindow()) {
            Rectangle visibleRect = context.getVisibleRect().intersection(rect);
            Graphics2D g2 = context.getGraphic2DForColor(Color.gray);
            String message = context.getReferenceFrame().getChrName().equals(Globals.CHR_ALL) ?
                    "Select a chromosome and zoom in to see alignments." :
                    "Zoom in to see alignments.";

            GraphicUtils.drawCenteredText(message, visibleRect, g2);
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
        downsampleRect.height = DOWNSAMPLED_ROW_HEIGHT;
        boolean downsampled = renderDownsampledIntervals(context, downsampleRect);

        alignmentsRect = new Rectangle(rect);
        if(downsampled) {
            alignmentsRect.y += DOWNSAMPLED_ROW_HEIGHT;
        }
        alignmentsRect.y +=  2;
        alignmentsRect.height -= (alignmentsRect.y - rect.y);
        renderAlignments(context, alignmentsRect);
    }

    private boolean renderDownsampledIntervals(RenderContext context, Rectangle downsampleRect) {

        // Might be offscreen
        if (!context.getVisibleRect().intersects(downsampleRect)) return false;

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
        if (renderOptions.getColorOption() == null && dataManager.hasYCTags()) {
            renderOptions.setColorOption(ColorOption.YC_TAG);
        }

        Map<String, PEStats> peStats = dataManager.getPEStats();
        if (peStats != null) {
            renderOptions.peStats = peStats;
        }

        Rectangle visibleRect = context.getVisibleRect();

        // Divide rectangle into equal height levels
        double y = inputRect.getY();
        double h;
        final DisplayMode displayMode = getDisplayMode();
        if (displayMode == DisplayMode.EXPANDED || displayMode == DisplayMode.FULL) {
            h = expandedHeight;
        } else if (displayMode == DisplayMode.COLLAPSED) {
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
                if (groupHeight > GROUP_LABEL_HEIGHT + 2 && !context.multiframe) {
                    String groupName = entry.getKey();
                    if (groupName.equals("SELECTED")) {
                        // Abbreviate the "SELECTED" group label to "S*" for display to conserve horizontal space in the UI.
                        groupName = "S*";
                    }
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

    /**
     * Render insertions at position of marker in expanded form, showing sequence.  Much of this method is just
     * a repeat of the loop for alignment rendering to compute the Y position of affected alignments.
     *
     * @param insertionMarker
     * @param context
     * @param inputRect
     */
    public void renderExpandedInsertion(InsertionMarker insertionMarker, RenderContext context, Rectangle inputRect) {

        boolean leaveMargin = getDisplayMode() != DisplayMode.SQUISHED;
        inputRect.y += DS_MARGIN_0 + DOWNSAMPLED_ROW_HEIGHT + DS_MARGIN_0;

        final AlignmentInterval loadedInterval = dataManager.getLoadedInterval(context.getReferenceFrame(), true);
        PackedAlignments groups = dataManager.getGroups(loadedInterval, renderOptions);
        if (groups == null) {
            //Assume we are still loading.
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
                    if (row.alignments != null)  // TODO -- not sure this is needed
                        BaseRenderer.drawExpandedInsertions(insertionMarker, row.alignments, context, rowRectangle, leaveMargin, renderOptions);
                    row.y = y;
                    row.h = h;
                }
                y += h;
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
            setDisplayMode(DisplayMode.SQUISHED);
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
            CircularViewUtilities.sendAlignmentsToJBrowse(inView, AlignmentTrack.this.getName(), chordColor);
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
            CircularViewUtilities.sendAlignmentsToJBrowse(inView, AlignmentTrack.this.getName(), chordColor);
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
        private static final Logger log = LogManager.getLogger(RenderOptions.class);

        private AlignmentTrack track;
        private Boolean shadeBasesOption;
        private Boolean shadeCenters;
        private Boolean flagUnmappedPairs;
        private Boolean showAllBases;
        private Integer minInsertSize;
        private Integer maxInsertSize;
        private ColorOption colorOption;
        private SortOption sortOption;
        private GroupOption groupByOption;
        private ShadeAlignmentsOption shadeAlignmentsOption;
        private DuplicatesOption duplicatesOption;
        private Integer mappingQualityLow;
        private Integer mappingQualityHigh;
        private boolean viewPairs = false;
        private String colorByTag;
        private String groupByTag;
        private String sortByTag;
        private String linkByTag;
        private Boolean linkedReads;
        private Boolean quickConsensusMode;
        private Boolean showMismatches;
        private Boolean indelQualColoring;
        private Boolean indelQualUsesMin;
        private Boolean indelQualSbx;
        private Boolean tailQualSbx;
        private Boolean hideTailSbx;
        private Boolean insertQualColoring;
        Boolean computeIsizes;
        private Double minInsertSizePercentile;
        private Double maxInsertSizePercentile;
        private Boolean pairedArcView;
        private Boolean flagZeroQualityAlignments;
        private Range groupByPos;
        private Boolean invertSorting;
        private boolean invertGroupSorting;
        private Boolean hideSmallIndels;
        private Integer smallIndelThreshold;
        private BaseModficationFilter basemodFilter;
        private Float basemodThreshold;
        private Boolean basemodDistinguishStrands;
        private int baseQualityMin;
        private int baseQualityMax;

        private Integer minJunctionCoverage;


        BisulfiteContext bisulfiteContext = BisulfiteContext.CG;
        Map<String, PEStats> peStats;

        RenderOptions(AlignmentTrack track) {
            this.track = track;
            peStats = new HashMap<>();

            // Set some constants -- for efficiency
            this.baseQualityMin = track == null ? 5 : track.getPreferences().getAsInt(SAM_BASE_QUALITY_MIN);
            this.baseQualityMax = track == null ? 20 : track.getPreferences().getAsInt(SAM_BASE_QUALITY_MAX);
        }

        IGVPreferences getPreferences() {
            return this.track != null ? this.track.getPreferences() : AlignmentTrack.getPreferences(ExperimentType.OTHER);
        }

        public int getMinJunctionCoverage() {
            return minJunctionCoverage != null ? minJunctionCoverage : PreferencesManager.getPreferences(Constants.RNA).getAsInt(SAM_JUNCTION_MIN_COVERAGE);
        }

        public void setMinJunctionCoverage(int minJunctionCoverage) {
            this.minJunctionCoverage = minJunctionCoverage;
        }

        public int getBaseQualityMin() {
            return baseQualityMin;
        }

        public int getBaseQualityMax() {
            return baseQualityMax;
        }

        public HashMap<String, Color> getSelectedReadNames() {
            return this.track.getSelectedReadNames();
        }

        public DisplayMode getDisplayMode() {
            return this.track.getDisplayMode();
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

        void setIndelQualColoring(boolean indelQualColoring) {
            this.indelQualColoring = indelQualColoring;
        }

        void setIndelQualUsesMin(boolean indelQualUsesMin) {
            this.indelQualUsesMin = indelQualUsesMin;
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

        void setSortOption(SortOption sortOption) {
            this.sortOption = sortOption;
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

        void setInvertSorting(boolean invertSorting) {
            this.invertSorting = invertSorting;
        }

        void setInvertGroupSorting(boolean invertGroupSorting) {
            this.invertGroupSorting = invertGroupSorting;
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

        void setShadeAlignmentsOption(ShadeAlignmentsOption shadeAlignmentsOption) {
            this.shadeAlignmentsOption = shadeAlignmentsOption;
        }

        void setShadeBasesOption(boolean shadeBasesOption) {
            this.shadeBasesOption = shadeBasesOption;
        }

        void setLinkedReads(boolean linkedReads) {
            this.linkedReads = linkedReads;
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

        public boolean isFlagZeroQualityAlignments() {
            return flagZeroQualityAlignments == null ? getPreferences().getAsBoolean(SAM_FLAG_ZERO_QUALITY) : flagZeroQualityAlignments;
        }

        public boolean isViewPairs() {
            return viewPairs;
        }

        public boolean isIndelQualColoring() {
            return indelQualColoring == null ? getPreferences().getAsBoolean(SAM_INDEL_QUAL_COLORING) : indelQualColoring;
        }

        public boolean isIndelQualUsesMin() {
            return indelQualUsesMin == null ? getPreferences().getAsBoolean(SAM_INDEL_QUAL_USES_MIN) : indelQualUsesMin;
        }


        // SBX Options
        public boolean isIndelQualSbx() {
            return ExperimentType.SBX == track.experimentType && (indelQualSbx == null ? getPreferences().getAsBoolean(SAM_INDEL_QUAL_SBX) : indelQualSbx);
        }

        public void setTailQualSbx(Boolean tailQualSbx) {
            this.tailQualSbx = tailQualSbx;
        }

        public boolean isTailQualSbx() {
            return ExperimentType.SBX == track.experimentType && (tailQualSbx == null ? getPreferences().getAsBoolean(SAM_TAIL_QUAL_SBX) : tailQualSbx);
        }

        public void setHideTailSbx(Boolean hideTailSbx) {
            this.hideTailSbx = hideTailSbx;
        }

        public boolean isHideTailSbx() {
            return ExperimentType.SBX == track.experimentType && (hideTailSbx == null ? getPreferences().getAsBoolean(SAM_HIDE_TAIL_SBX) : hideTailSbx);
        }

        public void setIndelQualSbx(Boolean indelQualSbx) {
            this.indelQualSbx = indelQualSbx;
        }
        // End SBX options

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

        public ShadeAlignmentsOption getShadeAlignmentsOption() {
            if (shadeAlignmentsOption != null) {
                return shadeAlignmentsOption;
            } else {
                try {
                    return ShadeAlignmentsOption.valueOf(getPreferences().get(SAM_SHADE_ALIGNMENT_BY));
                } catch (IllegalArgumentException e) {
                    log.error("Error parsing alignment shade option: " + ShadeAlignmentsOption.valueOf(getPreferences().get(SAM_SHADE_ALIGNMENT_BY)));
                    return ShadeAlignmentsOption.NONE;
                }
            }
        }

        public DuplicatesOption getDuplicatesOption() {
            final IGVPreferences prefs = getPreferences();
            if (duplicatesOption != null) {
                return duplicatesOption;
            } else {
                duplicatesOption = prefs.getAsBoolean(SAM_FILTER_DUPLICATES)
                        ? DuplicatesOption.FILTER
                        : DuplicatesOption.SHOW;
            }
            return duplicatesOption;
        }

        public void setDuplicatesOption(final DuplicatesOption duplicatesOption) {
            this.duplicatesOption = duplicatesOption;
        }

        public int getMappingQualityLow() {
            return mappingQualityLow == null ? getPreferences().getAsInt(SAM_SHADE_QUALITY_LOW) : mappingQualityLow;
        }

        public int getMappingQualityHigh() {
            return mappingQualityHigh == null ? getPreferences().getAsInt(SAM_SHADE_QUALITY_HIGH) : mappingQualityHigh;
        }

        SortOption getSortOption() {
            return sortOption == null ? SortOption.fromString(getPreferences().get(SAM_SORT_OPTION)) : sortOption;
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

        public boolean isInvertSorting() {
            return invertSorting == null ? getPreferences().getAsBoolean(SAM_INVERT_SORT) : invertSorting;
        }

        public boolean isInvertGroupSorting() {
            return invertGroupSorting;
        }

        public String getLinkByTag() {
            return linkByTag == null ? getPreferences().get(SAM_LINK_TAG) : linkByTag;
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
            return linkedReads == null ? getPreferences().getAsBoolean(SAM_LINKED_READS) : linkedReads;
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

        public BaseModficationFilter getBasemodFilter() {
            return basemodFilter;
        }

        public void setBasemodFilter(BaseModficationFilter basemodFilter) {

            this.basemodFilter = basemodFilter;
        }

        public float getBasemodThreshold() {
            return basemodThreshold == null ? getPreferences().getAsFloat(BASEMOD_THRESHOLD) : basemodThreshold.floatValue();
        }

        public void setBasemodThreshold(float basemodThreshold) {
            this.basemodThreshold = basemodThreshold;
        }

        public boolean getBasemodDistinguishStrands() {
            return basemodDistinguishStrands == null ? getPreferences().getAsBoolean(BASEMOD_DISTINGUISH_STRANDS) : basemodDistinguishStrands;
        }

        public void setBasemodDistinguishStrands(boolean basemodDistinguishStrands) {
            this.basemodDistinguishStrands = basemodDistinguishStrands;
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
            if (shadeAlignmentsOption != null) {
                element.setAttribute("shadeAlignmentsByOption", shadeAlignmentsOption.toString());
            }
            if (duplicatesOption != null) {
                element.setAttribute("duplicatesOption", duplicatesOption.toString());
            }
            if (mappingQualityLow != null) {
                element.setAttribute("mappingQualityLow", mappingQualityLow.toString());
            }
            if (mappingQualityHigh != null) {
                element.setAttribute("mappingQualityHigh", mappingQualityHigh.toString());
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
            if (invertSorting != null) {
                element.setAttribute("invertSorting", Boolean.toString(invertSorting));
            }
            if (sortOption != null) {
                element.setAttribute("sortOption", sortOption.toString());
            }
            if (invertGroupSorting) {
                element.setAttribute("invertGroupSorting", Boolean.toString(invertGroupSorting));
            }
            if (hideSmallIndels != null) {
                element.setAttribute("hideSmallIndels", hideSmallIndels.toString());
            }
            if (smallIndelThreshold != null) {
                element.setAttribute("smallIndelThreshold", smallIndelThreshold.toString());
            }
            if (basemodFilter != null) {
                element.setAttribute("basemodFilter", basemodFilter.toString());
            }
            if (basemodThreshold != null) {
                element.setAttribute("basemodThredhold", String.valueOf(basemodThreshold));
            }
            if (minJunctionCoverage != null) {
                element.setAttribute("minJunctionCoverage", String.valueOf(minJunctionCoverage));
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
                // Convert deprecated options
                final String attributeValue = element.getAttribute("colorOption");
                if ("BASE_MODIFICATION_6MA".equals(attributeValue)) {
                    colorOption = ColorOption.BASE_MODIFICATION;
                    basemodFilter = new BaseModficationFilter("a");
                } else if ("BASE_MODIFICATION_5MC".equals(attributeValue)) {
                    colorOption = ColorOption.BASE_MODIFICATION_2COLOR;
                    // basemodFilter = new BaseModficationFilter(null, 'C');
                } else if ("BASE_MODIFICATION_C".equals(attributeValue)) {
                    colorOption = ColorOption.BASE_MODIFICATION;
                    // basemodFilter = new BaseModficationFilter(null, 'C');
                } else {
                    colorOption = ColorOption.valueOf(attributeValue);
                }
            }
            if (element.hasAttribute("sortOption")) {
                sortOption = SortOption.fromString((element.getAttribute("sortOption")));
            }
            if (element.hasAttribute("groupByOption")) {
                String value = element.getAttribute("groupByOption");
                if (value.equals("HAPLOTYPE")) {
                    value = "CLUSTER";  // Backward compatibility
                }
                groupByOption = GroupOption.valueOf(value);
            }
            if (element.hasAttribute("shadeAlignmentsByOption")) {
                shadeAlignmentsOption = ShadeAlignmentsOption.valueOf(element.getAttribute("shadeAlignmentsByOption"));
            }
            if (element.hasAttribute("duplicatesOption")) {
                duplicatesOption = CollUtils.valueOf(DuplicatesOption.class, element.getAttribute("duplicatesOption"), null);
            }
            if (element.hasAttribute("mappingQualityLow")) {
                mappingQualityLow = Integer.parseInt(element.getAttribute("mappingQualityLow"));
            }
            if (element.hasAttribute("mappingQualityHigh")) {
                mappingQualityHigh = Integer.parseInt(element.getAttribute("mappingQualityHigh"));
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
            if (element.hasAttribute("invertSorting")) {
                invertSorting = Boolean.parseBoolean(element.getAttribute("invertSorting"));
            }
            if (element.hasAttribute("invertGroupSorting")) {
                invertGroupSorting = Boolean.parseBoolean(element.getAttribute("invertGroupSorting"));
            }
            if (element.hasAttribute("hideSmallIndels")) {
                hideSmallIndels = Boolean.parseBoolean(element.getAttribute("hideSmallIndels"));
            }
            if (element.hasAttribute("smallIndelThreshold")) {
                smallIndelThreshold = Integer.parseInt(element.getAttribute("smallIndelThreshold"));
            }
            if (element.hasAttribute("showInsertionMarkers")) {
                // TODO -- something with this
                // showInsertionMarkers = Boolean.parseBoolean(element.getAttribute("showInsertionMarkers"));
            }
            if (element.hasAttribute("basemodFilter")) {
                basemodFilter = BaseModficationFilter.fromString(element.getAttribute("basemodFilter"));
            }
            if (element.hasAttribute("basemodThreshold")) {
                basemodFilter = BaseModficationFilter.fromString(element.getAttribute("basemodThreshold"));
            }
            if (element.hasAttribute("minJunctionCoverage")) {
                minJunctionCoverage = Integer.parseInt(element.getAttribute("minJunctionCoverage"));
            }
        }

    }


}
