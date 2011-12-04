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
import org.broad.igv.ui.event.TrackGroupEvent;
import org.broad.igv.ui.panel.*;
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
import java.lang.ref.SoftReference;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AlignmentTrack extends AbstractTrack implements AlignmentTrackEventListener {


    public enum SortOption {
        START, STRAND, NUCELOTIDE, QUALITY, SAMPLE, READ_GROUP, INSERT_SIZE
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


    public enum ColorOption {
        INSERT_SIZE, READ_STRAND, FRAGMENT_STRAND, PAIR_ORIENTATION, SAMPLE, READ_GROUP, BISULFITE, NOMESEQ, NONE;
    }


    public static final int MIN_ALIGNMENT_SPACING = 10;
    static final ColorOption DEFAULT_COLOR_OPTION = ColorOption.INSERT_SIZE;
    static final boolean DEFAULT_SHOWALLBASES = false;
    static final BisulfiteContext DEFAULT_BISULFITE_CONTEXT = BisulfiteContext.CG;

    private static ColorOption colorByOption = null;
    private static BisulfiteContext bisulfiteContext = null;

    private SequenceTrack sequenceTrack;
    private CoverageTrack coverageTrack;
    private SpliceJunctionFinderTrack spliceJunctionTrack;

    private RenderOptions renderOptions;

    private static Logger log = Logger.getLogger(AlignmentTrack.class);
    private int expandedHeight = 14;
    private int maxCollapsedHeight = 4;
    private int collapsedHeight = maxCollapsedHeight;

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

        renderer = new AlignmentRenderer();

        this.setDisplayMode(DisplayMode.EXPANDED);

        if (prefs.getAsBoolean(PreferenceManager.SAM_SHOW_REF_SEQ)) {
            sequenceTrack = new SequenceTrack("Reference sequence");
            sequenceTrack.setHeight(14);
        }

        renderOptions = new RenderOptions();


        if (colorByOption == null) {
            String colorByString = PreferenceManager.getInstance().get(PreferenceManager.SAM_COLOR_BY);
            if (colorByString == null) {
                colorByOption = DEFAULT_COLOR_OPTION;
            } else {
                try {
                    colorByOption = ColorOption.valueOf(colorByString);
                } catch (Exception e) {
                    log.error("Error setting color option", e);
                    colorByOption = DEFAULT_COLOR_OPTION;
                }
            }
        }
        // Override color option preference, if necessary
        if (!dataManager.isPairedEnd() &&
                (colorByOption == ColorOption.INSERT_SIZE || colorByOption == ColorOption.PAIR_ORIENTATION)) {
            colorByOption = ColorOption.NONE;
        }

        if (bisulfiteContext == null) {
            String str = PreferenceManager.getInstance().get(PreferenceManager.SAM_BISULFITE_CONTEXT);
            if (str == null) {
                bisulfiteContext = DEFAULT_BISULFITE_CONTEXT;
            } else {
                try {
                    bisulfiteContext = BisulfiteContext.valueOf(str);
                } catch (Exception e) {
                    log.error("Error setting bisulfite option", e);
                    bisulfiteContext = DEFAULT_BISULFITE_CONTEXT;

                }
            }
        }

        // Register track
        IGV.getInstance().addAlignmentTrackEventListener(this);


    }


    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public void setSpliceJunctionTrack(SpliceJunctionFinderTrack spliceJunctionTrack) {
        this.spliceJunctionTrack = spliceJunctionTrack;
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
        int h = Math.max(minHeight, getNLevels() * (getRowHeight()) + 20);
        return h;
    }

    private int getRowHeight() {
        return getDisplayMode() == DisplayMode.EXPANDED ? expandedHeight : collapsedHeight;
    }

    private int getNLevels() {
        return dataManager.getNLevels();
    }

    public void render(RenderContext context, Rectangle rect) {

        parent = context.getPanel();

        // Split rects
        int seqHeight = sequenceTrack == null ? 0 : sequenceTrack.getHeight();
        if (seqHeight > 0) {
            Rectangle seqRect = new Rectangle(rect);
            seqRect.height = seqHeight;
            sequenceTrack.render(context, seqRect);
        }

        int gap = (seqHeight > 0 ? seqHeight : 6);

        rect.y += gap;
        rect.height -= gap;
        renderedRect = new Rectangle(rect);

        if (context.getScale() > minVisibleScale) {

            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Zoom in to see alignments.", context.getVisibleRect(), g);
            return;

        }

        renderFeatures(context, rect);
    }

    private void renderFeatures(RenderContext context, Rectangle inputRect) {
        try {
            log.debug("Render features");
            List<AlignmentInterval.Row> tmp = dataManager.getAlignmentRows(context); //genomeId, chr, start, end);

            Map<String, PEStats> peStats = dataManager.getPEStats();
            if (peStats != null) {
                renderOptions.peStats = peStats;
            }

            if (tmp == null) {
                return;
            }


            Rectangle visibleRect = context.getVisibleRect();

            // Divide rectangle into equal height levels
            double y = inputRect.getY();
            double h = expandedHeight;
            if (getDisplayMode() != DisplayMode.EXPANDED) {
                int visHeight = context.getVisibleRect().height;
                int depth = dataManager.getMaxLevels();
                collapsedHeight = Math.min(maxCollapsedHeight, Math.max(1, Math.min(expandedHeight, visHeight / depth)));
                h = collapsedHeight;
            }

            int levelNumber = 0;
            for (AlignmentInterval.Row row : tmp) {

                if ((visibleRect != null && y > visibleRect.getMaxY())) {
                    return;
                }

                if (y + h > visibleRect.getY()) {
                    Rectangle rect = new Rectangle(inputRect.x, (int) y, inputRect.width, (int) h);
                    renderOptions.colorOption = colorByOption;
                    renderOptions.bisulfiteContextRenderOption = bisulfiteContext;
                    renderer.renderAlignments(row.alignments,
                            context,
                            rect,
                            renderOptions,
                            getDisplayMode() == DisplayMode.EXPANDED,
                            selectedReadNames);
                }
                levelNumber++;
                y += h;
            }
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

    public void sortRows(SortOption option, ReferenceFrame referenceFrame, double location) {
        dataManager.sortRows(option, referenceFrame, location);
    }

    /**
     * Sort alignment rows such that alignments that intersect from the
     * center appear left to right by start position
     */
    public void groupRows(SortOption option, ReferenceFrame referenceFrame) {
        dataManager.repackAlignments(referenceFrame, option);
    }


    public void packAlignments(ReferenceFrame referenceFrame) {
        dataManager.repackAlignments(referenceFrame);
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

        Alignment feature = getAlignmentAt(position, y, frame);
        if (feature == null) {
            return null;
        }
        return feature.getValueString(position, getWindowFunction());
    }

    private Alignment getAlignmentAt(double position, int y, ReferenceFrame frame) {

        List<AlignmentInterval.Row> alignmentRows = dataManager.getAlignmentRows(frame);

        if (alignmentRows == null || alignmentRows.isEmpty()) {
            return null;
        }

        int h = getDisplayMode() == DisplayMode.EXPANDED ? expandedHeight : collapsedHeight;
        int levelNumber = (y - renderedRect.y) / h;
        if (levelNumber < 0 || levelNumber >= alignmentRows.size()) {
            return null;
        }

        AlignmentInterval.Row row = alignmentRows.get(levelNumber);
        List<Alignment> features = row.alignments;

        // No buffer for alignments,  you must zoom in far enough for them to be visible
        int buffer = 0;
        return (Alignment) FeatureUtils.getFeatureAt(position, buffer, features);

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


    private static Alignment getFeatureContaining(
            List<Alignment> features, int right) {

        int leftBounds = 0;
        int rightBounds = features.size() - 1;
        int idx = features.size() / 2;
        int lastIdx = -1;

        while (idx != lastIdx) {
            lastIdx = idx;
            Alignment f = features.get(idx);
            if (f.contains(right)) {
                return f;
            }

            if (f.getStart() > right) {
                rightBounds = idx;
                idx = (leftBounds + idx) / 2;
            } else {
                leftBounds = idx;
                idx = (rightBounds + idx) / 2;

            }

        }
        // Check the extremes
        if (features.get(0).contains(right)) {
            return features.get(0);
        }

        if (features.get(rightBounds).contains(right)) {
            return features.get(rightBounds);
        }

        return null;
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
    boolean showCenterLine;
    boolean flagUnmappedPairs;
    boolean showAllBases;
    private boolean computeIsizes;
    private int minInsertSize;
    private int maxInsertSize;
    private double minInsertSizePercentile;
    private double maxInsertSizePercentile;
    ColorOption colorOption;
    BisulfiteContext bisulfiteContextRenderOption;
    //ContinuousColorScale insertSizeColorScale;
    private boolean viewPairs = false;
    public boolean flagZeroQualityAlignments = true;
    Map<String, PEStats> peStats;

    RenderOptions() {
        PreferenceManager prefs = PreferenceManager.getInstance();
        shadeBases = prefs.getAsBoolean(PreferenceManager.SAM_SHADE_BASE_QUALITY);
        shadeCenters = prefs.getAsBoolean(PreferenceManager.SAM_SHADE_CENTER);
        showCenterLine = prefs.getAsBoolean(PreferenceManager.SAM_SHOW_CENTER_LINE);
        flagUnmappedPairs = prefs.getAsBoolean(PreferenceManager.SAM_FLAG_UNMAPPED_PAIR);
        computeIsizes = prefs.getAsBoolean(PreferenceManager.SAM_COMPUTE_ISIZES);
        minInsertSize = prefs.getAsInt(PreferenceManager.SAM_MIN_INSERT_SIZE_THRESHOLD);
        maxInsertSize = prefs.getAsInt(PreferenceManager.SAM_MAX_INSERT_SIZE_THRESHOLD);
        minInsertSizePercentile = prefs.getAsFloat(PreferenceManager.SAM_MIN_INSERT_SIZE_PERCENTILE);
        maxInsertSizePercentile = prefs.getAsFloat(PreferenceManager.SAM_MAX_INSERT_SIZE_PERCENTILE);
        showAllBases = DEFAULT_SHOWALLBASES;
        colorOption = colorByOption;
        bisulfiteContextRenderOption = bisulfiteContext;
        flagZeroQualityAlignments = prefs.getAsBoolean(PreferenceManager.SAM_FLAG_ZERO_QUALITY);
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
        if (showCenterLine != prefs.getAsBoolean(PreferenceManager.SAM_SHOW_CENTER_LINE)) {
            attributes.put("shadeCenters", String.valueOf(showCenterLine));
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
            attributes.put("colorOption", colorByOption.toString());
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
            colorByOption = colorOption;
        }
        value = attributes.get("bisulfiteContextRenderOption");
        if (value != null) {
            bisulfiteContextRenderOption = BisulfiteContext.valueOf(value);
            bisulfiteContext = bisulfiteContextRenderOption;
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
        addPackMenuItem();
        addCoverageDepthMenuItem();

        addSeparator();
        addColorByMenuItem();
        addShadeBaseMenuItem();
        addShowAllBasesMenuItem();

        addSeparator();
        addViewAsPairsMenuItem();
        addGoToMate(e);
        showMateRegion(e);
        addInsertSizeMenuItem();

        addSeparator();
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
            nomeESeqOption.setSelected(colorByOption == ColorOption.NOMESEQ);
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
            m1.setSelected(bisulfiteContext == item);
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

        JMenuItem m2 = new JMenuItem("by strand");
        m2.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {
                IGV.getInstance().sortAlignmentTracks(SortOption.STRAND);
                refresh();

            }
        });
        groupMenu.add(m2);

        JMenuItem m5 = new JMenuItem("by sample");
        m5.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {

                IGV.getInstance().sortAlignmentTracks(SortOption.SAMPLE);
                refresh();

            }
        });
        groupMenu.add(m5);

        JMenuItem m6 = new JMenuItem("by read group");
        m6.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {

                IGV.getInstance().sortAlignmentTracks(SortOption.READ_GROUP);
                refresh();

            }
        });
        groupMenu.add(m6);

        add(groupMenu);
    }

    public void addSortMenuItem() {//ReferenceFrame frame) {
        // Change track height by attribute
        JMenu item = new JMenu("Sort alignments");

        JMenuItem m1 = new JMenuItem("by start location");
        m1.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {
                IGV.getInstance().sortAlignmentTracks(SortOption.START);
                refresh();

            }
        });
        item.add(m1);

        JMenuItem m2 = new JMenuItem("by strand");
        m2.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {
                IGV.getInstance().sortAlignmentTracks(SortOption.STRAND);
                refresh();

            }
        });
        item.add(m2);

        JMenuItem m3 = new JMenuItem("by base");
        m3.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {

                IGV.getInstance().sortAlignmentTracks(SortOption.NUCELOTIDE);
                refresh();

            }
        });
        item.add(m3);

        JMenuItem m4 = new JMenuItem("by mapping quality");
        m4.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {

                IGV.getInstance().sortAlignmentTracks(SortOption.QUALITY);
                refresh();

            }
        });
        item.add(m4);


        JMenuItem m5 = new JMenuItem("by sample");
        m5.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {

                IGV.getInstance().sortAlignmentTracks(SortOption.SAMPLE);
                refresh();

            }
        });
        item.add(m5);

        JMenuItem m6 = new JMenuItem("by read group");
        m6.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent aEvt) {

                IGV.getInstance().sortAlignmentTracks(SortOption.READ_GROUP);
                refresh();

            }
        });
        item.add(m6);

        if (dataManager.isPairedEnd()) {
            JMenuItem m7 = new JMenuItem("by insert size");
            m7.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {

                    IGV.getInstance().sortAlignmentTracks(SortOption.INSERT_SIZE);
                    refresh();

                }
            });
            item.add(m7);
        }


        add(item);
    }


    private void setBisulfiteContext(BisulfiteContext option) {
        bisulfiteContext = option;
        PreferenceManager.getInstance().put(PreferenceManager.SAM_BISULFITE_CONTEXT, option.toString());
    }

    private void setColorOption(ColorOption option) {
        colorByOption = option;
        PreferenceManager.getInstance().put(PreferenceManager.SAM_COLOR_BY, option.toString());
    }

    public void addColorByMenuItem() {
        // Change track height by attribute
        JMenu colorMenu = new JMenu("Color alignments");

        ButtonGroup group = new ButtonGroup();

        JRadioButtonMenuItem noneOption = new JRadioButtonMenuItem("no color");
        noneOption.setSelected(colorByOption == ColorOption.NONE);
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
            isizeOption.setSelected(colorByOption == ColorOption.INSERT_SIZE);
            isizeOption.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent aEvt) {
                    setColorOption(ColorOption.INSERT_SIZE);
                    refresh();
                }
            });
            colorMenu.add(isizeOption);
            group.add(isizeOption);

            JRadioButtonMenuItem m1a = new JRadioButtonMenuItem("by pair orientation");
            m1a.setSelected(colorByOption == ColorOption.PAIR_ORIENTATION);
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
        readStrandOption.setSelected(colorByOption == ColorOption.READ_STRAND);
        readStrandOption.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aEvt) {
                setColorOption(ColorOption.READ_STRAND);
                refresh();
            }
        });
        colorMenu.add(readStrandOption);
        group.add(readStrandOption);

        JRadioButtonMenuItem readGroupOption = new JRadioButtonMenuItem("by read group");
        readGroupOption.setSelected(colorByOption == ColorOption.READ_GROUP);
        readGroupOption.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aEvt) {
                setColorOption(ColorOption.READ_GROUP);
                refresh();
            }
        });
        colorMenu.add(readGroupOption);
        group.add(readGroupOption);

        JRadioButtonMenuItem sampleOption = new JRadioButtonMenuItem("by sample");
        sampleOption.setSelected(colorByOption == ColorOption.SAMPLE);
        sampleOption.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aEvt) {
                setColorOption(ColorOption.SAMPLE);
                refresh();
            }
        });
        colorMenu.add(sampleOption);
        group.add(sampleOption);

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
        item.setSelected(dataManager.isLoadAsPairs());
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent aEvt) {

                dataManager.setLoadAsPairs(item.isSelected());
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

        if (colorByOption == ColorOption.BISULFITE || colorByOption == ColorOption.NOMESEQ) {
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
