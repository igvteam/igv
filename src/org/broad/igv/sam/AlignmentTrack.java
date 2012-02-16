/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
package org.broad.igv.sam;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.goby.GobyCountArchiveDataSource;
import org.broad.igv.lists.GeneList;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.session.Session;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.InsertSizeSettingsDialog;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.event.AlignmentTrackEvent;
import org.broad.igv.ui.event.AlignmentTrackEventListener;
import org.broad.igv.ui.panel.DataPanel;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.Pair;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AlignmentTrack extends AbstractTrack implements AlignmentTrackEventListener {

    private static Logger log = Logger.getLogger(AlignmentTrack.class);
    static final int GROUP_MARGIN = 5;
    static final int TOP_MARGIN = 20;


    public enum ColorOption {
        INSERT_SIZE, READ_STRAND, FIRST_OF_PAIR_STRAND, PAIR_ORIENTATION, SAMPLE, READ_GROUP, BISULFITE, NOMESEQ, TAG, NONE;
    }

    public enum SortOption {
        START, STRAND, NUCELOTIDE, QUALITY, SAMPLE, READ_GROUP, INSERT_SIZE, FIRST_OF_PAIR_STRAND, MATE_CHR, TAG
    }

    public enum GroupOption {
        STRAND, SAMPLE, READ_GROUP, FIRST_OF_PAIR_STRAND, TAG, NONE
    }

    public enum BisulfiteContext {
        CG, CHH, CHG, HCG, GCH, WCG
    }

    protected static final Map<BisulfiteContext, String> bisulfiteContextToPubString = new HashMap<BisulfiteContext, String>();

    static {
        bisulfiteContextToPubString.put(BisulfiteContext.CG, "CG");
        bisulfiteContextToPubString.put(BisulfiteContext.CHH, "CHH");
        bisulfiteContextToPubString.put(BisulfiteContext.CHG, "CHG");
        bisulfiteContextToPubString.put(BisulfiteContext.HCG, "HCG");
        bisulfiteContextToPubString.put(BisulfiteContext.GCH, "GCH");
        bisulfiteContextToPubString.put(BisulfiteContext.WCG, "WCG");
    }

    protected static final Map<BisulfiteContext, Pair<byte[], byte[]>> bisulfiteContextToContextString = new HashMap<BisulfiteContext, Pair<byte[], byte[]>>();

    static {
        bisulfiteContextToContextString.put(BisulfiteContext.CG, new Pair<byte[], byte[]>(new byte[]{}, new byte[]{'G'}));
        bisulfiteContextToContextString.put(BisulfiteContext.CHH, new Pair<byte[], byte[]>(new byte[]{}, new byte[]{'H', 'H'}));
        bisulfiteContextToContextString.put(BisulfiteContext.CHG, new Pair<byte[], byte[]>(new byte[]{}, new byte[]{'H', 'G'}));
        bisulfiteContextToContextString.put(BisulfiteContext.HCG, new Pair<byte[], byte[]>(new byte[]{'H'}, new byte[]{'G'}));
        bisulfiteContextToContextString.put(BisulfiteContext.GCH, new Pair<byte[], byte[]>(new byte[]{'G'}, new byte[]{'H'}));
        bisulfiteContextToContextString.put(BisulfiteContext.WCG, new Pair<byte[], byte[]>(new byte[]{'W'}, new byte[]{'G'}));
    }

    public static final int MIN_ALIGNMENT_SPACING = 10;
    static final ColorOption DEFAULT_COLOR_OPTION = ColorOption.INSERT_SIZE;
    static final boolean DEFAULT_SHOWALLBASES = false;
    static final BisulfiteContext DEFAULT_BISULFITE_CONTEXT = BisulfiteContext.CG;

    //private static GroupOption groupByOption = null;
    //private static ColorOption colorByOption = null;
    //private static BisulfiteContext bisulfiteContext = null;

    private SequenceTrack sequenceTrack;
    private CoverageTrack coverageTrack;
    private SpliceJunctionFinderTrack spliceJunctionTrack;

    private RenderOptions renderOptions;

    private int expandedHeight = 14;
    private int maxSquishedHeight = 4;
    private int squishedHeight = maxSquishedHeight;

    private FeatureRenderer renderer;
    private double minVisibleScale = 25;
    private Rectangle renderedRect;
    private HashMap<String, Color> selectedReadNames = new HashMap();
    private int selectionColorIndex = 0;
    private int minHeight = 100;
    private AlignmentDataManager dataManager;
    Genome genome;

    // The "parent" of the track (a DataPanel).  This release of IGV does not support owner-track releationships
    // directory,  so this field might be null at any given time.  It is updated each repaint.
    DataPanel parent;


    /**
     * Create a new alignment track
     *
     * @param locator
     * @param dataManager
     * @param genome
     */
    public AlignmentTrack(ResourceLocator locator, AlignmentDataManager dataManager, Genome genome) {
        super(locator);

        this.genome = genome;
        this.dataManager = dataManager;

        PreferenceManager prefs = PreferenceManager.getInstance();

        float maxRange = prefs.getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        minVisibleScale = (maxRange * 1000) / 700;

        renderer = AlignmentRenderer.getInstance();

        this.setDisplayMode(DisplayMode.EXPANDED);

        if (prefs.getAsBoolean(PreferenceManager.SAM_SHOW_REF_SEQ)) {
            sequenceTrack = new SequenceTrack("Reference sequence");
            sequenceTrack.setHeight(14);
        }

        renderOptions = new RenderOptions();

        // Register track
        if (!Globals.isHeadless()) {
            IGV.getInstance().addAlignmentTrackEventListener(this);
        }

    }


    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public void setSpliceJunctionTrack(SpliceJunctionFinderTrack spliceJunctionTrack) {
        this.spliceJunctionTrack = spliceJunctionTrack;
        dataManager.setSpliceJunctionTrack(spliceJunctionTrack);
    }


    public SpliceJunctionFinderTrack getSpliceJunctionTrack() {
        return spliceJunctionTrack;
    }

    public void setRenderer(FeatureRenderer renderer) {
        this.renderer = renderer;
    }

    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {
        return new PopupMenu(te);
    }

    @Override
    public void setHeight(int preferredHeight) {
        super.setHeight(preferredHeight);
        minHeight = preferredHeight;
    }

    @Override
    public int getHeight() {
        // TODO -- what is the 20 for? JTR
        int nGroups = dataManager.getMaxGroupCount();
        int h = Math.max(minHeight, getNLevels() * getRowHeight() + nGroups * GROUP_MARGIN +
                TOP_MARGIN);
        return h;
    }

    private int getRowHeight() {
        return getDisplayMode() == DisplayMode.EXPANDED ? expandedHeight : squishedHeight;
    }

    private int getNLevels() {
        return dataManager.getNLevels();
    }

    public void render(RenderContext context, Rectangle rect) {

        parent = context.getPanel();

        // Split track rectangle into sections.
        int seqHeight = sequenceTrack == null ? 0 : sequenceTrack.getHeight();
        if (seqHeight > 0) {
            Rectangle seqRect = new Rectangle(rect);
            seqRect.height = seqHeight;
            sequenceTrack.render(context, seqRect);
        }

        // Top gap.  If there's a sequence track no gap is needed
        int gap = (seqHeight > 0 ? seqHeight : 6);

        rect.y += gap;
        rect.height -= gap;
        renderedRect = new Rectangle(rect);

        if (context.getScale() > minVisibleScale) {
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Zoom in to see alignments.", context.getVisibleRect(), g);
            return;

        }

        renderAlignments(context, rect);
    }

    private void renderAlignments(RenderContext context, Rectangle inputRect) {
        try {
            log.debug("Render features");
            Map<String, List<AlignmentInterval.Row>> groups =
                    dataManager.getGroups(context, renderOptions.groupByOption, renderOptions.getGroupByTag());

            Map<String, PEStats> peStats = dataManager.getPEStats();
            if (peStats != null) {
                renderOptions.peStats = peStats;
            }

            if (groups == null) {
                return;
            }


            Rectangle visibleRect = context.getVisibleRect();
            final boolean leaveMargin = getDisplayMode() == DisplayMode.EXPANDED;

            // Divide rectangle into equal height levels
            double y = inputRect.getY();
            double h = expandedHeight;
            if (getDisplayMode() != DisplayMode.EXPANDED) {
                int visHeight = context.getVisibleRect().height;
                int depth = dataManager.getMaxLevels();
                squishedHeight = Math.min(maxSquishedHeight, Math.max(1, Math.min(expandedHeight, visHeight / depth)));
                h = squishedHeight;
            }

            // Loop through groups
            Graphics2D groupBorderGraphics = context.getGraphic2DForColor(AlignmentRenderer.GROUP_DIVIDER_COLOR);
            int nGroups = groups.size();
            int groupNumber = 0;
            for (Map.Entry<String, List<AlignmentInterval.Row>> entry : groups.entrySet()) {
                String group = entry.getKey();
                groupNumber++;

                // Loop through the alignment rows for this group
                List<AlignmentInterval.Row> rows = entry.getValue();
                for (AlignmentInterval.Row row : rows) {

                    if ((visibleRect != null && y > visibleRect.getMaxY())) {
                        return;
                    }

                    if (y + h > visibleRect.getY()) {
                        Rectangle rowRectangle = new Rectangle(inputRect.x, (int) y, inputRect.width, (int) h);
                        renderer.renderAlignments(row.alignments,
                                context,
                                rowRectangle,
                                renderOptions,
                                leaveMargin,
                                selectedReadNames);
                    }
                    y += h;
                }

                // Draw a subtle divider line between groups
                if (groupNumber < nGroups) {
                    int borderY = (int) y + GROUP_MARGIN / 2;
                    groupBorderGraphics.drawLine(inputRect.x, borderY, inputRect.width, borderY);
                }
                y += GROUP_MARGIN;
            }

            final int bottom = inputRect.y + inputRect.height;
            groupBorderGraphics.drawLine(inputRect.x, bottom, inputRect.width, bottom);

        } catch (Exception ex) {
            log.error("Error rendering track", ex);
            throw new RuntimeException("Error rendering track ", ex);

        }

    }

    public void clearCaches() {
        dataManager.clear();
        renderOptions = new RenderOptions();
    }


    /**
     * Sort alignment rows based on alignments that intersent location
     */

    public void sortRows(SortOption option, ReferenceFrame referenceFrame, double location, String tag) {
        dataManager.sortRows(option, referenceFrame, location, tag);
    }

    /**
     * Sort alignment rows such that alignments that intersect from the
     * center appear left to right by start position
     */
    public void groupAlignments(GroupOption option, ReferenceFrame referenceFrame) {
        if (renderOptions.groupByOption != option) {
            renderOptions.groupByOption = (option == GroupOption.NONE ? null : option);
            dataManager.repackAlignments(referenceFrame, option, renderOptions.getGroupByTag());
        }
    }


    public void packAlignments(ReferenceFrame referenceFrame) {
        dataManager.repackAlignments(referenceFrame, renderOptions.groupByOption, renderOptions.getGroupByTag());
    }

    /**
     * Copy the contents of the popup text to the system clipboard.
     */
    public void copyToClipboard(final TrackClickEvent e, Alignment alignment, double location) {

        if (alignment != null) {
            StringBuffer buf = new StringBuffer();
            buf.append(alignment.getValueString(location, null).replace("<br>", "\n"));
            buf.append("\n");
            buf.append("Alignment start position = " + alignment.getChr() + ":" + (alignment.getAlignmentStart() + 1));
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
    public void gotoMate(final TrackClickEvent te, Alignment alignment) {


        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {

                setSelected(alignment);

                String chr = mate.getChr();
                int start = mate.start - 1;
                te.getFrame().centerOnLocation(chr, start);
                te.getFrame().recordHistory();
            } else {
                MessageUtils.showMessage("Alignment does not have mate, or it is not mapped.");
            }
        }
    }

    /**
     * Split the screen so the current view and mate region are side by side.  Need a better
     * name for this method.
     */
    public void splitScreenMate(final TrackClickEvent te, Alignment alignment) {

        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {

                setSelected(alignment);

                String mateChr = mate.getChr();
                int mateStart = mate.start - 1;

                ReferenceFrame frame = te.getFrame();
                String locus1 = frame.getFormattedLocusString();

                // Generate a locus string for the read mate.  Keep the window width (in base pairs) == to the current range
                ReferenceFrame.Range range = frame.getCurrentRange();
                int length = range.getLength();
                int s2 = Math.max(0, mateStart - length / 2);
                int e2 = s2 + length;
                String startStr = NumberFormat.getInstance().format(s2);
                String endStr = NumberFormat.getInstance().format(e2);
                String mateLocus = mateChr + ":" + startStr + "-" + endStr;

                Session currentSession = IGV.getInstance().getSession();

                // If we are already in gene list mode add the mate as another panel, otherwise switch to gl mode

                List<String> loci = null;
                if (FrameManager.isGeneListMode()) {
                    loci = new ArrayList(FrameManager.getFrames().size());
                    for (ReferenceFrame ref : FrameManager.getFrames()) {
                        loci.add(ref.getInitialLocus().toString());
                    }
                    loci.add(mateLocus);
                } else {
                    loci = Arrays.asList(locus1, mateLocus);
                }
                GeneList.sortByPosition(loci);
                StringBuffer listName = new StringBuffer();
                for (String s : loci) {
                    listName.append(s + "   ");
                }

                GeneList geneList = new GeneList(listName.toString(), loci, false);
                currentSession.setCurrentGeneList(geneList);
                IGV.getInstance().resetFrames();

            } else {
                MessageUtils.showMessage("Alignment does not have mate, or it is not mapped.");
            }
        }
    }


    public void setWindowFunction(WindowFunction type) {
        // ignored
    }

    public WindowFunction getWindowFunction() {
        return null;
    }

    public void setRendererClass(Class rc) {
        // ignored
    }

    // SamTracks use a custom renderer, not derived from Renderer

    public Renderer getRenderer() {
        return null;
    }

    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, ReferenceFrame frame) {
        return 0.0f;
    }

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        Alignment feature = getAlignmentAt(position + 1, y, frame);
        if (feature == null) {
            return null;
        }
        return feature.getValueString(position, getWindowFunction());
    }

    private Alignment getAlignmentAt(double position, int y, ReferenceFrame frame) {

        Map<String, List<AlignmentInterval.Row>> groups = dataManager.getGroupedAlignments(frame);

        if (groups == null || groups.isEmpty()) {
            return null;
        }

        int h = getDisplayMode() == DisplayMode.EXPANDED ? expandedHeight : squishedHeight;
//        int levelNumber = (y - renderedRect.y) / h;
//        if (levelNumber < 0) {
//            return null;
//        }

        int startY = renderedRect.y;
        final boolean leaveMargin = getDisplayMode() == DisplayMode.EXPANDED;

        for (List<AlignmentInterval.Row> rows : groups.values()) {
            int endY = startY + rows.size() * h;
            if (y >= startY && y < endY) {
                int levelNumber = (y - startY) / h;
                AlignmentInterval.Row row = rows.get(levelNumber);
                List<Alignment> features = row.alignments;

                // No buffer for alignments,  you must zoom in far enough for them to be visible
                int buffer = 0;
                return (Alignment) FeatureUtils.getFeatureAt(position, buffer, features);
            }
            startY = endY + GROUP_MARGIN;
        }

        return null;

    }

    /**
     * Handle an AlignmentTrackEvent.
     *
     * @param e
     */
    public void onAlignmentTrackEvent(AlignmentTrackEvent e) {

        AlignmentTrackEvent.Type type = e.getType();
        switch (type) {
            case VISIBILITY_WINDOW:
                visibilityWindowChanged();
                break;
            case RELOAD:
            case SPLICE_JUNCTION:
                clearCaches();
                break;
        }

    }

    /**
     * The visibility window has changed.
     */
    private void visibilityWindowChanged() {
        PreferenceManager prefs = PreferenceManager.getInstance();
        float maxRange = prefs.getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        minVisibleScale = (maxRange * 1000) / 700;
    }


    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();
        if (Globals.IS_MAC && e.isMetaDown() || (!Globals.IS_MAC && e.isControlDown())) {
            // Selection
            final ReferenceFrame frame = te.getFrame();
            if (frame != null) {
                selectAlignment(e, frame);
                if (parent != null) parent.repaint();
                return true;
            }

        }
        return false;
    }

    private void selectAlignment(MouseEvent e, ReferenceFrame frame) {
        double location = frame.getChromosomePosition(e.getX());
        double displayLocation = location + 1;
        Alignment alignment = this.getAlignmentAt(displayLocation, e.getY(), frame);
        if (alignment != null) {
            if (selectedReadNames.containsKey(alignment.getReadName())) {
                selectedReadNames.remove(alignment.getReadName());
            } else {
                setSelected(alignment);
            }

        }

    }

    private void setSelected(Alignment alignment) {
        Color c = alignment.isPaired() && alignment.getMate() != null && alignment.getMate().isMapped() ?
                ColorUtilities.randomColor(selectionColorIndex++) : Color.black;
        selectedReadNames.put(alignment.getReadName(), c);
    }


    private void refresh() {
        IGV.getInstance().getContentPane().getMainPanel().invalidate();
        IGV.getInstance().repaintDataPanels();
    }

    public static boolean isBisulfiteColorType(ColorOption o) {
        return (o.equals(ColorOption.BISULFITE) || o.equals(ColorOption.NOMESEQ));
    }

    public static String getBisulfiteContextPubStr(BisulfiteContext item) {
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

    @Override
    public Map<String, String> getPersistentState() {
        Map<String, String> attrs = super.getPersistentState();
        attrs.putAll(renderOptions.getPersistentState());
        return attrs;
    }

    @Override
    public void restorePersistentState(Map<String, String> attributes) {
        super.restorePersistentState(attributes);    //To change body of overridden methods use File | Settings | File Templates.
        renderOptions.restorePersistentState(attributes);
    }


    public static class RenderOptions {
        boolean shadeBases;
        boolean shadeCenters;
        boolean flagUnmappedPairs;
        boolean showAllBases;
        boolean showMismatches = true;
        private boolean computeIsizes;
        private int minInsertSize;
        private int maxInsertSize;
        private double minInsertSizePercentile;
        private double maxInsertSizePercentile;
        private ColorOption colorOption;
        GroupOption groupByOption = null;
        BisulfiteContext bisulfiteContext;
        //ContinuousColorScale insertSizeColorScale;
        private boolean viewPairs = false;
        public boolean flagZeroQualityAlignments = true;

        Map<String, PEStats> peStats;

        private String colorByTag;
        private String groupByTag;
        private String sortByTag;

        RenderOptions() {
            PreferenceManager prefs = PreferenceManager.getInstance();
            shadeBases = prefs.getAsBoolean(PreferenceManager.SAM_SHADE_BASE_QUALITY);
            shadeCenters = prefs.getAsBoolean(PreferenceManager.SAM_SHADE_CENTER);
            flagUnmappedPairs = prefs.getAsBoolean(PreferenceManager.SAM_FLAG_UNMAPPED_PAIR);
            computeIsizes = prefs.getAsBoolean(PreferenceManager.SAM_COMPUTE_ISIZES);
            minInsertSize = prefs.getAsInt(PreferenceManager.SAM_MIN_INSERT_SIZE_THRESHOLD);
            maxInsertSize = prefs.getAsInt(PreferenceManager.SAM_MAX_INSERT_SIZE_THRESHOLD);
            minInsertSizePercentile = prefs.getAsFloat(PreferenceManager.SAM_MIN_INSERT_SIZE_PERCENTILE);
            maxInsertSizePercentile = prefs.getAsFloat(PreferenceManager.SAM_MAX_INSERT_SIZE_PERCENTILE);
            showAllBases = DEFAULT_SHOWALLBASES;
            colorOption = ColorOption.valueOf(prefs.get(PreferenceManager.SAM_COLOR_BY));
            GroupOption groupByOption = null;
            flagZeroQualityAlignments = prefs.getAsBoolean(PreferenceManager.SAM_FLAG_ZERO_QUALITY);
            bisulfiteContext = DEFAULT_BISULFITE_CONTEXT;

            colorByTag = prefs.get(PreferenceManager.SAM_COLOR_BY_TAG);
            sortByTag = prefs.get(PreferenceManager.SAM_SORT_BY_TAG);
            groupByTag = prefs.get(PreferenceManager.SAM_GROUP_BY_TAG);

            //updateColorScale();
        }

        /*private void updateColorScale() {
            int delta = 1;
            if (medianInsertSize == 0 || madInsertSize == 0) {
                delta = (maxInsertSizeThreshold - minInsertSizeThreshold) / 10;
            } else {
                delta = Math.min((maxInsertSizeThreshold - minInsertSizeThreshold) / 3, madInsertSize);
            }
            insertSizeColorScale = new ContinuousColorScale(minInsertSizeThreshold, minInsertSizeThreshold + delta,
                    maxInsertSizeThreshold - delta, maxInsertSizeThreshold, Color.blue, AlignmentRenderer.grey1, Color.red);
        }*/

        /**
         * Called by session writer.  Return instance variable values as a map of strings.  Used to record current state
         * of object.   Variables with default values are not stored, as it is presumed the user has not changed them.
         *
         * @return
         */
        public Map<String, String> getPersistentState() {
            Map<String, String> attributes = new HashMap();
            PreferenceManager prefs = PreferenceManager.getInstance();
            if (shadeBases != prefs.getAsBoolean(PreferenceManager.SAM_SHADE_BASE_QUALITY)) {
                attributes.put("shadeBases", String.valueOf(shadeBases));
            }
            if (shadeCenters != prefs.getAsBoolean(PreferenceManager.SAM_SHADE_CENTER)) {
                attributes.put("shadeCenters", String.valueOf(shadeBases));
            }
            if (flagUnmappedPairs != prefs.getAsBoolean(PreferenceManager.SAM_FLAG_UNMAPPED_PAIR)) {
                attributes.put("flagUnmappedPairs", String.valueOf(flagUnmappedPairs));
            }
            if (maxInsertSize != prefs.getAsInt(PreferenceManager.SAM_MAX_INSERT_SIZE_THRESHOLD)) {
                attributes.put("insertSizeThreshold", String.valueOf(maxInsertSize));
            }
            if (getMinInsertSize() != prefs.getAsInt(PreferenceManager.SAM_MIN_INSERT_SIZE_THRESHOLD)) {
                attributes.put("minInsertSizeThreshold", String.valueOf(maxInsertSize));
            }
            if (showAllBases != DEFAULT_SHOWALLBASES) {
                attributes.put("showAllBases", String.valueOf(showAllBases));
            }
            if (colorOption != DEFAULT_COLOR_OPTION) {
                attributes.put("colorOption", colorOption.toString());
            }
            if (groupByOption != null) {
                attributes.put("groupByOption", groupByOption.toString());
            }
            if (colorByTag != null && colorByTag.length() > 0) {
                attributes.put("colorByTag", colorByTag);
            }
            if (groupByTag != null && groupByTag.length() > 0) {
                attributes.put("groupByTag", groupByTag);
            }
            if (sortByTag != null) {
                attributes.put("sortByTag", sortByTag);
            }

            return attributes;
        }

        /**
         * Called by session reader.  Restores state of object.
         *
         * @param attributes
         */
        public void restorePersistentState(Map<String, String> attributes) {

            String value;
            value = attributes.get("insertSizeThreshold");
            if (value != null) {
                maxInsertSize = Integer.parseInt(value);
            }
            value = attributes.get("minInsertSizeThreshold");
            if (value != null) {
                setMinInsertSize(Integer.parseInt(value));
            }
            value = attributes.get("shadeBases");
            if (value != null) {
                shadeBases = Boolean.parseBoolean(value);
            }
            value = attributes.get("shadeCenters");
            if (value != null) {
                shadeCenters = Boolean.parseBoolean(value);
            }
            value = attributes.get("flagUnmappedPairs");
            if (value != null) {
                flagUnmappedPairs = Boolean.parseBoolean(value);
            }
            value = attributes.get("showAllBases");
            if (value != null) {
                showAllBases = Boolean.parseBoolean(value);
            }
            value = attributes.get("colorOption");
            if (value != null) {
                colorOption = ColorOption.valueOf(value);
            }
            value = attributes.get("groupByOption");
            if (value != null) {
                groupByOption = GroupOption.valueOf(value);
            }
            value = attributes.get("bisulfiteContextRenderOption");
            if (value != null) {
                bisulfiteContext = BisulfiteContext.valueOf(value);
            }
            value = attributes.get("groupByTag");
            if (value != null) {
                groupByTag = value;
            }
            value = attributes.get("sortByTag");
            if (value != null) {
                sortByTag = value;
            }
            value = attributes.get("colorByTag");
            if (value != null) {
                colorByTag = value;
            }
        }

        public int getMinInsertSize() {
            return minInsertSize;
        }

        public void setMinInsertSize(int minInsertSize) {
            this.minInsertSize = minInsertSize;
            //updateColorScale();
        }

        public int getMaxInsertSize() {
            return maxInsertSize;

        }


        public boolean isViewPairs() {
            return viewPairs;
        }

        public void setViewPairs(boolean viewPairs) {
            this.viewPairs = viewPairs;
        }

        public boolean isComputeIsizes() {
            return computeIsizes;
        }

        public void setComputeIsizes(boolean computeIsizes) {
            this.computeIsizes = computeIsizes;
        }

        public double getMinInsertSizePercentile() {
            return minInsertSizePercentile;
        }

        public void setMinInsertSizePercentile(double minInsertSizePercentile) {
            this.minInsertSizePercentile = minInsertSizePercentile;
        }

        public double getMaxInsertSizePercentile() {
            return maxInsertSizePercentile;
        }

        public void setMaxInsertSizePercentile(double maxInsertSizePercentile) {
            this.maxInsertSizePercentile = maxInsertSizePercentile;
        }

        public void setMaxInsertSize(int maxInsertSize) {
            this.maxInsertSize = maxInsertSize;
        }

        public ColorOption getColorOption() {
            return colorOption;
        }

        public void setColorOption(ColorOption colorOption) {
            this.colorOption = colorOption;
        }

        public void setColorByTag(String colorByTag) {
            this.colorByTag = colorByTag;
            PreferenceManager.getInstance().put(PreferenceManager.SAM_COLOR_BY_TAG, colorByTag);
        }

        public String getColorByTag() {
            return colorByTag;
        }

        public String getSortByTag() {
            return sortByTag;
        }

        public void setSortByTag(String sortByTag) {
            this.sortByTag = sortByTag;
        }

        public String getGroupByTag() {
            return groupByTag;
        }

        public void setGroupByTag(String groupByTag) {
            this.groupByTag = groupByTag;
        }
    }

    class PopupMenu extends IGVPopupMenu {

        PopupMenu(final TrackClickEvent e) {

            Collection<Track> tracks = new ArrayList();
            tracks.add(AlignmentTrack.this);

            JLabel popupTitle = new JLabel("  " + AlignmentTrack.this.getName(), JLabel.CENTER);

            Font newFont = getFont().deriveFont(Font.BOLD, 12);
            popupTitle.setFont(newFont);
            if (popupTitle != null) {
                add(popupTitle);
            }
            addSeparator();
            add(TrackMenuUtils.getTrackRenameItem(tracks));
            addCopyToClipboardItem(e);

            addSeparator();
            addGroupMenuItem();
            addSortMenuItem();
            addColorByMenuItem();

            addSeparator();
            addShadeBaseMenuItem();
            addShowMismatchesMenuItem();
            addShowAllBasesMenuItem();

            addSeparator();
            addViewAsPairsMenuItem();
            addGoToMate(e);
            showMateRegion(e);
            addInsertSizeMenuItem();

            addSeparator();
            addPackMenuItem();
            addCoverageDepthMenuItem();
            addShowCoverageItem();
            addLoadCoverageDataItem();

            addSeparator();
            TrackMenuUtils.addDisplayModeItems(tracks, this);

            addSeparator();
            addSelecteByNameItem();
            addClearSelectionsMenuItem();

            addSeparator();
            add(TrackMenuUtils.getRemoveMenuItem(tracks));

            return;
        }

        private JMenu getBisulfiteContextMenuItem(ButtonGroup group) {
            // Change track height by attribute
            //JMenu bisulfiteContextMenu = new JMenu("Bisulfite Contexts");
            JMenu bisulfiteContextMenu = new JMenu("bisulfite mode");


            JRadioButtonMenuItem nomeESeqOption = null;
            boolean showNomeESeq = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_NOMESEQ_ENABLED);
            if (showNomeESeq) {
                nomeESeqOption = new JRadioButtonMenuItem("NOMe-seq bisulfite mode");
                nomeESeqOption.setSelected(renderOptions.colorOption == ColorOption.NOMESEQ);
                nomeESeqOption.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent aEvt) {
                        setColorOption(ColorOption.NOMESEQ);
                        refresh();
                    }
                });
                group.add(nomeESeqOption);
            }

            for (final BisulfiteContext item : BisulfiteContext.values()) {

                String optionStr = getBisulfiteContextPubStr(item);
                JRadioButtonMenuItem m1 = new JRadioButtonMenuItem(optionStr);
                m1.setSelected(renderOptions.bisulfiteContext == item);
                m1.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent aEvt) {
                        setColorOption(ColorOption.BISULFITE);
                        setBisulfiteContext(item);
                        refresh();
                    }
                });
                bisulfiteContextMenu.add(m1);
                group.add(m1);
            }

            if (nomeESeqOption != null) {
                bisulfiteContextMenu.add(nomeESeqOption);
            }

            return bisulfiteContextMenu;

        }

        public void addSelecteByNameItem() {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("Select by name...");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    String val = MessageUtils.showInputDialog("Enter read name: ");
                    if (val != null && val.trim().length() > 0) {
                        selectedReadNames.put(val, ColorUtilities.randomColor(selectedReadNames.size() + 1));
                        refresh();
                    }
                }
            });

            add(item);
        }

        public void addGroupMenuItem() {//ReferenceFrame frame) {
            // Change track height by attribute
            JMenu groupMenu = new JMenu("Group alignments");
            ButtonGroup group = new ButtonGroup();

            JCheckBoxMenuItem m1 = new JCheckBoxMenuItem("none");
            m1.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    IGV.getInstance().groupAlignmentTracks(GroupOption.NONE);
                    refresh();

                }
            });
            m1.setSelected(renderOptions.groupByOption == null);
            groupMenu.add(m1);
            group.add(m1);

            JCheckBoxMenuItem m2 = new JCheckBoxMenuItem("by read strand");
            m2.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    IGV.getInstance().groupAlignmentTracks(GroupOption.STRAND);
                    refresh();

                }
            });
            m2.setSelected(renderOptions.groupByOption == GroupOption.STRAND);
            groupMenu.add(m2);
            group.add(m2);

            JCheckBoxMenuItem fragmentStrandOption = new JCheckBoxMenuItem("by first-of-pair strand");
            fragmentStrandOption.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    IGV.getInstance().groupAlignmentTracks(GroupOption.FIRST_OF_PAIR_STRAND);
                    refresh();

                }
            });
            fragmentStrandOption.setSelected(renderOptions.groupByOption == GroupOption.FIRST_OF_PAIR_STRAND);
            groupMenu.add(fragmentStrandOption);
            group.add(fragmentStrandOption);

            JCheckBoxMenuItem m5 = new JCheckBoxMenuItem("by sample");
            m5.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {

                    IGV.getInstance().groupAlignmentTracks(GroupOption.SAMPLE);
                    refresh();

                }
            });
            m5.setSelected(renderOptions.groupByOption == GroupOption.SAMPLE);
            groupMenu.add(m5);
            group.add(m5);

            JCheckBoxMenuItem m6 = new JCheckBoxMenuItem("by read group");
            m6.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {

                    IGV.getInstance().groupAlignmentTracks(GroupOption.READ_GROUP);
                    refresh();

                }
            });
            m6.setSelected(renderOptions.groupByOption == GroupOption.READ_GROUP);
            groupMenu.add(m6);
            group.add(m6);


            JCheckBoxMenuItem tagOption = new JCheckBoxMenuItem("by tag");
            tagOption.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getGroupByTag());
                    renderOptions.setGroupByTag(tag);
                    IGV.getInstance().groupAlignmentTracks(GroupOption.TAG);
                    refresh();

                }
            });
            tagOption.setSelected(renderOptions.groupByOption == GroupOption.TAG);
            groupMenu.add(tagOption);
            group.add(tagOption);


//            JMenuItem tagOption = new JMenuItem("by tag");
//            tagOption.addActionListener(new ActionListener() {
//                public void actionPerformed(ActionEvent aEvt) {
//                    String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getTagKey());
//                    renderOptions.setTagKey(tag);
//                    IGV.getInstance().sortAlignmentTracks(SortOption.TAG, tag);
//                    refresh();
//                }
//            });
//            sortMenu.add(tagOption);
//

            add(groupMenu);
        }


        /**
         * Sort menu
         */
        public void addSortMenuItem() {


            JMenu sortMenu = new JMenu("Sort alignments");
            JMenuItem m1 = new JMenuItem("by start location");
            m1.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    IGV.getInstance().sortAlignmentTracks(SortOption.START, null);
                    refresh();

                }
            });
            sortMenu.add(m1);

            JMenuItem m2 = new JMenuItem("by read strand");
            m2.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    IGV.getInstance().sortAlignmentTracks(SortOption.STRAND, null);
                    refresh();

                }
            });
            sortMenu.add(m2);

            JMenuItem fragmentStrandOption = new JMenuItem("by first-of-pair strand");
            fragmentStrandOption.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    IGV.getInstance().sortAlignmentTracks(SortOption.FIRST_OF_PAIR_STRAND, null);
                    refresh();

                }
            });
            sortMenu.add(fragmentStrandOption);

            JMenuItem m3 = new JMenuItem("by base");
            m3.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {

                    IGV.getInstance().sortAlignmentTracks(SortOption.NUCELOTIDE, null);
                    refresh();

                }
            });
            sortMenu.add(m3);

            JMenuItem m4 = new JMenuItem("by mapping quality");
            m4.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {

                    IGV.getInstance().sortAlignmentTracks(SortOption.QUALITY, null);
                    refresh();

                }
            });
            sortMenu.add(m4);


            JMenuItem m5 = new JMenuItem("by sample");
            m5.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {

                    IGV.getInstance().sortAlignmentTracks(SortOption.SAMPLE, null);
                    refresh();

                }
            });
            sortMenu.add(m5);

            JMenuItem m6 = new JMenuItem("by read group");
            m6.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {

                    IGV.getInstance().sortAlignmentTracks(SortOption.READ_GROUP, null);
                    refresh();

                }
            });
            sortMenu.add(m6);

            if (dataManager.isPairedEnd()) {
                JMenuItem m7 = new JMenuItem("by insert size");
                m7.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {

                        IGV.getInstance().sortAlignmentTracks(SortOption.INSERT_SIZE, null);
                        refresh();

                    }
                });
                sortMenu.add(m7);
            }

            if (dataManager.isPairedEnd()) {
                JMenuItem m7 = new JMenuItem("by chromosome of mate");
                m7.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {

                        IGV.getInstance().sortAlignmentTracks(SortOption.MATE_CHR, null);
                        refresh();

                    }
                });
                sortMenu.add(m7);
            }

            JMenuItem tagOption = new JMenuItem("by tag");
            tagOption.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getSortByTag());
                    renderOptions.setSortByTag(tag);
                    IGV.getInstance().sortAlignmentTracks(SortOption.TAG, tag);
                    refresh();
                }
            });
            sortMenu.add(tagOption);


            add(sortMenu);
        }


        private void setBisulfiteContext(BisulfiteContext option) {
            renderOptions.bisulfiteContext = option;
            PreferenceManager.getInstance().put(PreferenceManager.SAM_BISULFITE_CONTEXT, option.toString());
        }

        private void setColorOption(ColorOption option) {
            renderOptions.colorOption = option;
            PreferenceManager.getInstance().put(PreferenceManager.SAM_COLOR_BY, option.toString());
        }

        public void addColorByMenuItem() {
            // Change track height by attribute
            JMenu colorMenu = new JMenu("Color alignments");

            ButtonGroup group = new ButtonGroup();

            JRadioButtonMenuItem noneOption = new JRadioButtonMenuItem("no color");
            noneOption.setSelected(renderOptions.colorOption == ColorOption.NONE);
            noneOption.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    setColorOption(ColorOption.NONE);
                    refresh();
                }
            });
            colorMenu.add(noneOption);
            group.add(noneOption);


            if (dataManager.isPairedEnd()) {
                JRadioButtonMenuItem isizeOption = new JRadioButtonMenuItem("by insert size");
                isizeOption.setSelected(renderOptions.colorOption == ColorOption.INSERT_SIZE);
                isizeOption.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent aEvt) {
                        setColorOption(ColorOption.INSERT_SIZE);
                        refresh();
                    }
                });
                colorMenu.add(isizeOption);
                group.add(isizeOption);

                JRadioButtonMenuItem m1a = new JRadioButtonMenuItem("by pair orientation");
                m1a.setSelected(renderOptions.colorOption == ColorOption.PAIR_ORIENTATION);
                m1a.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent aEvt) {
                        setColorOption(ColorOption.PAIR_ORIENTATION);
                        refresh();
                    }
                });
                colorMenu.add(m1a);
                group.add(m1a);
            }

            JRadioButtonMenuItem readStrandOption = new JRadioButtonMenuItem("by read strand");
            readStrandOption.setSelected(renderOptions.colorOption == ColorOption.READ_STRAND);
            readStrandOption.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    setColorOption(ColorOption.READ_STRAND);
                    refresh();
                }
            });
            colorMenu.add(readStrandOption);
            group.add(readStrandOption);

            if (dataManager.isPairedEnd()) {
                JRadioButtonMenuItem fragmentStrandOption = new JRadioButtonMenuItem("by first-of-pair strand");
                fragmentStrandOption.setSelected(renderOptions.colorOption == ColorOption.FIRST_OF_PAIR_STRAND);
                fragmentStrandOption.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent aEvt) {
                        setColorOption(ColorOption.FIRST_OF_PAIR_STRAND);
                        refresh();
                    }
                });
                colorMenu.add(fragmentStrandOption);
                group.add(fragmentStrandOption);
            }

            JRadioButtonMenuItem readGroupOption = new JRadioButtonMenuItem("by read group");
            readGroupOption.setSelected(renderOptions.colorOption == ColorOption.READ_GROUP);
            readGroupOption.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    setColorOption(ColorOption.READ_GROUP);
                    refresh();
                }
            });
            colorMenu.add(readGroupOption);
            group.add(readGroupOption);

            JRadioButtonMenuItem sampleOption = new JRadioButtonMenuItem("by sample");
            sampleOption.setSelected(renderOptions.colorOption == ColorOption.SAMPLE);
            sampleOption.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    setColorOption(ColorOption.SAMPLE);
                    refresh();
                }
            });
            colorMenu.add(sampleOption);
            group.add(sampleOption);


            JRadioButtonMenuItem tagOption = new JRadioButtonMenuItem("by tag");
            tagOption.setSelected(renderOptions.colorOption == ColorOption.TAG);
            tagOption.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    setColorOption(ColorOption.TAG);
                    String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getColorByTag());
                    renderOptions.setColorByTag(tag);
                    PreferenceManager.getInstance();
                    refresh();
                }
            });
            colorMenu.add(tagOption);
            group.add(tagOption);


            colorMenu.add(getBisulfiteContextMenuItem(group));


            add(colorMenu);

        }


        public void addPackMenuItem() {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("Re-pack alignments");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    UIUtilities.invokeOnEventThread(new Runnable() {

                        public void run() {
                            IGV.getInstance().packAlignmentTracks();
                            refresh();
                        }
                    });
                }
            });

            add(item);
        }

        public void addCopyToClipboardItem(final TrackClickEvent te) {

            final MouseEvent me = te.getMouseEvent();
            JMenuItem item = new JMenuItem("Copy read details to clipboard");

            final ReferenceFrame frame = te.getFrame();
            if (frame == null) {
                item.setEnabled(false);
            } else {
                final double location = frame.getChromosomePosition(me.getX());
                double displayLocation = location + 1;
                final Alignment alignment = getAlignmentAt(displayLocation, me.getY(), frame);

                // Change track height by attribute
                item.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {
                        copyToClipboard(te, alignment, location);

                    }
                });
                if (alignment == null) {
                    item.setEnabled(false);
                }
            }

            add(item);
        }

        public void addViewAsPairsMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("View as pairs");
            item.setSelected(dataManager.isViewAsPairs());
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    boolean viewAsPairs = item.isSelected();

                    // TODO -- generalize this test to all incompatible pairings
                    if (viewAsPairs && renderOptions.groupByOption == GroupOption.STRAND) {
                        boolean ungroup = MessageUtils.confirm("\"View as pairs\" is incompatible with \"Group by strand\". Ungroup?");
                        if (ungroup) {
                            renderOptions.groupByOption = null;
                        } else {
                            return;
                        }
                    }

                    dataManager.setViewAsPairs(item.isSelected(), renderOptions.groupByOption, renderOptions.getGroupByTag());
                    refresh();
                }
            });
            item.setEnabled(dataManager.isPairedEnd());
            add(item);
        }

        public void addGoToMate(final TrackClickEvent te) {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("Go to mate");
            MouseEvent e = te.getMouseEvent();

            final ReferenceFrame frame = te.getFrame();
            if (frame == null) {
                item.setEnabled(false);
            } else {
                double location = frame.getChromosomePosition(e.getX());

                double displayLocation = location + 1;
                final Alignment alignment = getAlignmentAt(displayLocation, e.getY(), frame);
                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent aEvt) {
                        gotoMate(te, alignment);
                    }
                });
                if (alignment == null || !alignment.isPaired() || !alignment.getMate().isMapped()) {
                    item.setEnabled(false);
                }
            }
            add(item);
        }


        public void showMateRegion(final TrackClickEvent te) {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("View mate region in split screen");
            MouseEvent e = te.getMouseEvent();

            final ReferenceFrame frame = te.getFrame();
            if (frame == null) {
                item.setEnabled(false);
            } else {
                double location = frame.getChromosomePosition(e.getX());
                double displayLocation = location + 1;
                final Alignment alignment = getAlignmentAt(displayLocation, e.getY(), frame);

                item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent aEvt) {
                        splitScreenMate(te, alignment);
                    }
                });
                if (alignment == null || !alignment.isPaired() || !alignment.getMate().isMapped()) {
                    item.setEnabled(false);
                }
            }
            add(item);
        }

        public void addClearSelectionsMenuItem() {
            // Change track height by attribute
            JMenuItem item = new JMenuItem("Clear selections");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    selectedReadNames.clear();
                    refresh();
                }
            });
            add(item);
        }

        public void addShowAllBasesMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Show all bases");

            if (renderOptions.colorOption == ColorOption.BISULFITE || renderOptions.colorOption == ColorOption.NOMESEQ) {
                //    item.setEnabled(false);
            } else {
                item.setSelected(renderOptions.showAllBases);
            }
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    renderOptions.showAllBases = item.isSelected();
                    refresh();
                }
            });
            add(item);
        }

        public void addShowMismatchesMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Show mismatched bases");


            item.setSelected(renderOptions.showMismatches);
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    renderOptions.showMismatches = item.isSelected();
                    refresh();
                }
            });
            add(item);
        }

        public void addCoverageDepthMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Set maximum coverage depth ...");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    int maxLevels = dataManager.getMaxLevels();
                    String val = MessageUtils.showInputDialog("Maximum coverage depth", String.valueOf(maxLevels));
                    try {
                        int newMaxLevels = Integer.parseInt(val);
                        if (newMaxLevels != maxLevels) {
                            dataManager.setMaxLevels(newMaxLevels);
                            //dataManager.reload();
                            refresh();
                        }
                    } catch (NumberFormatException ex) {
                        MessageUtils.showMessage("Insert size must be an integer value: " + val);
                    }

                }
            });
            add(item);
        }

        public void addInsertSizeMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Set insert size options ...");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {

                    InsertSizeSettingsDialog dlg = new InsertSizeSettingsDialog(IGV.getMainFrame(), renderOptions);
                    dlg.setModal(true);
                    dlg.setVisible(true);
                    if (!dlg.isCanceled()) {
                        renderOptions.setComputeIsizes(dlg.isComputeIsize());
                        renderOptions.setMinInsertSizePercentile(dlg.getMinPercentile());
                        renderOptions.setMaxInsertSizePercentile(dlg.getMaxPercentile());
                        if (renderOptions.isComputeIsizes()) {
                            dataManager.updatePEStats(renderOptions);
                        }

                        renderOptions.setMinInsertSize(dlg.getMinThreshold());
                        renderOptions.setMaxInsertSize(dlg.getMaxThreshold());
                        refresh();
                    }
                }
            });


            item.setEnabled(dataManager.isPairedEnd());
            add(item);
        }


        public void addShadeBaseMenuItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Shade base by quality");
            item.setSelected(renderOptions.shadeBases);
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    UIUtilities.invokeOnEventThread(new Runnable() {

                        public void run() {
                            renderOptions.shadeBases = item.isSelected();
                            refresh();
                        }
                    });
                }
            });

            add(item);
        }


        public void addShowCoverageItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Show coverage track");
            item.setSelected(getCoverageTrack() != null && getCoverageTrack().isVisible());
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    UIUtilities.invokeOnEventThread(new Runnable() {

                        public void run() {
                            if (getCoverageTrack() != null) {
                                getCoverageTrack().setVisible(item.isSelected());
                                refresh();
                                IGV.getInstance().repaintNamePanels();
                            }
                        }
                    });
                }
            });

            add(item);
        }

        public void addLoadCoverageDataItem() {
            // Change track height by attribute
            final JMenuItem item = new JCheckBoxMenuItem("Load coverage data...");
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    UIUtilities.invokeOnEventThread(new Runnable() {

                        public void run() {

                            final PreferenceManager prefs = PreferenceManager.getInstance();
                            File initDirectory = prefs.getLastTrackDirectory();
                            File file = FileDialogUtils.chooseFile("Select coverage file", initDirectory, FileDialog.LOAD);
                            if (file != null) {
                                prefs.setLastTrackDirectory(file.getParentFile());
                                String path = file.getAbsolutePath();
                                if (path.endsWith(".tdf") || path.endsWith(".tdf")) {
                                    TDFReader reader = TDFReader.getReader(file.getAbsolutePath());
                                    TDFDataSource ds = new TDFDataSource(reader, 0, getName() + " coverage", genome);
                                    getCoverageTrack().setDataSource(ds);
                                    refresh();
                                } else if (path.endsWith(".counts")) {
                                    DataSource ds = new GobyCountArchiveDataSource(file);
                                    getCoverageTrack().setDataSource(ds);
                                    refresh();
                                } else {
                                    MessageUtils.showMessage("Coverage data must be in .tdf format");
                                }
                            }
                        }
                    });
                }
            });

            add(item);
        }

    }


}
