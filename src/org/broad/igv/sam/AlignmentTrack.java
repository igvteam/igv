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
package org.broad.igv.sam;

//~--- non-JDK imports --------------------------------------------------------

import com.iontorrent.data.FlowDistribution;
import com.iontorrent.data.ReadInfo;
import com.iontorrent.utils.LocationListener;
import com.iontorrent.utils.SimpleDialog;
import com.iontorrent.views.FlowSignalDistributionPanel;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.data.CoverageDataSource;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.goby.GobyCountArchiveDataSource;
import org.broad.igv.lists.GeneList;
import org.broad.igv.renderer.GraphicUtils;
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
import org.broad.igv.util.collections.CollUtils;

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
    static final int DS_MARGIN_0 = 2;
    static final int DOWNAMPLED_ROW_HEIGHT = 3;
    static final int DS_MARGIN_2 = 5;

    public enum ShadeBasesOption {
        NONE, QUALITY, FLOW_SIGNAL_DEVIATION_READ, FLOW_SIGNAL_DEVIATION_REFERENCE
    }

    public enum ExperimentType {
        RNA, BISULFITE, OTHER
    }

    public enum ColorOption {
        INSERT_SIZE, READ_STRAND, FIRST_OF_PAIR_STRAND, PAIR_ORIENTATION, SAMPLE, READ_GROUP, BISULFITE, NOMESEQ,
        TAG, NONE, UNEXPECTED_PAIR
    }

    public enum SortOption {
        START, STRAND, NUCELOTIDE, QUALITY, SAMPLE, READ_GROUP, INSERT_SIZE, FIRST_OF_PAIR_STRAND, MATE_CHR, TAG;
    }

    public enum GroupOption {
        STRAND, SAMPLE, READ_GROUP, FIRST_OF_PAIR_STRAND, TAG, PAIR_ORIENTATION, MATE_CHROMOSOME, NONE
    }

    public enum BisulfiteContext {
        CG, CHH, CHG, HCG, GCH, WCG
    }

    enum OrientationType {
        RR, LL, RL, LR, UNKNOWN
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

    static final ShadeBasesOption DEFAULT_SHADE_BASES_OPTION = ShadeBasesOption.QUALITY;
    static final ColorOption DEFAULT_COLOR_OPTION = ColorOption.INSERT_SIZE;
    static final boolean DEFAULT_SHOWALLBASES = false;
    static final BisulfiteContext DEFAULT_BISULFITE_CONTEXT = BisulfiteContext.CG;
    private boolean ionTorrent;
    private SequenceTrack sequenceTrack;
    private CoverageTrack coverageTrack;
    private SpliceJunctionFinderTrack spliceJunctionTrack;
    private final RenderOptions renderOptions = new RenderOptions();
    private int expandedHeight = 14;
    private int maxSquishedHeight = 4;
    private int squishedHeight = maxSquishedHeight;
    private FeatureRenderer renderer;
    private double minVisibleScale = 25;
    private HashMap<String, Color> selectedReadNames = new HashMap();
    private int selectionColorIndex = 0;
    private int minHeight = 50;
    private AlignmentDataManager dataManager;
    private Genome genome;
    // The "parent" of the track (a DataPanel).  This field might be null at any given time.  It is updated each repaint.
    JComponent parent;
    private Rectangle alignmentsRect;
    private Rectangle downsampleRect;

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

        ionTorrent = dataManager.isIonTorrent();

        minimumHeight = 50;
        maximumHeight = Integer.MAX_VALUE;

        PreferenceManager prefs = PreferenceManager.getInstance();

        dataManager.setShowSpliceJunctions(prefs.getAsBoolean(PreferenceManager.SAM_SHOW_JUNCTION_TRACK));

        float maxRange = prefs.getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        minVisibleScale = (maxRange * 1000) / 700;

        renderer = AlignmentRenderer.getInstance();

        this.setDisplayMode(DisplayMode.EXPANDED);

        if (prefs.getAsBoolean(PreferenceManager.SAM_SHOW_REF_SEQ)) {
            sequenceTrack = new SequenceTrack("Reference sequence");
            sequenceTrack.setHeight(14);
        }

        if (renderOptions.getColorOption() == ColorOption.BISULFITE) {
            setExperimentType(ExperimentType.BISULFITE);
        }

        // Register track
        if (!Globals.isHeadless()) {
            IGV.getInstance().addAlignmentTrackEventListener(this);
        }

    }

    @Override
    public void updateGenome(Genome genome) {
        dataManager.updateGenome(genome);
    }

    /**
     * Set the experiment type (RNA, Bisulfite, or OTHER)
     *
     * @param type
     */
    public void setExperimentType(ExperimentType type) {
        dataManager.setExperimentType(type);
        if (spliceJunctionTrack != null) {
            spliceJunctionTrack.setVisible(type != ExperimentType.BISULFITE);
        }
    }

    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
        this.coverageTrack.setRenderOptions(this.renderOptions);
    }

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }

    public void setSpliceJunctionTrack(SpliceJunctionFinderTrack spliceJunctionTrack) {
        this.spliceJunctionTrack = spliceJunctionTrack;
        if (dataManager.getExperimentType() == ExperimentType.BISULFITE) {
            spliceJunctionTrack.setVisible(false);
        }
    }


    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {
        return new PopupMenu(te);
    }

    @Override
    public void setHeight(int preferredHeight) {
        super.setHeight(preferredHeight);
        minimumHeight = preferredHeight;
    }

    @Override
    public int getHeight() {

        if (parent != null
                && (parent instanceof DataPanel && ((DataPanel) parent).getFrame().getScale() > minVisibleScale)) {
            return minimumHeight;
        }

        int nGroups = dataManager.getMaxGroupCount();

        int h = Math.max(minHeight, getNLevels() * getRowHeight() + nGroups * GROUP_MARGIN + TOP_MARGIN
                + DS_MARGIN_0 + DOWNAMPLED_ROW_HEIGHT + DS_MARGIN_2);


        h = Math.min(maximumHeight, h);
        return h;
    }

    private int getRowHeight() {
        return getDisplayMode() == DisplayMode.EXPANDED ? expandedHeight : squishedHeight;
    }

    private int getNLevels() {
        return dataManager.getNLevels();
    }

    @Override
    public void preload(RenderContext context) {
        dataManager.preload(context, renderOptions, true);
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

        // Top gap.
        rect.y += DS_MARGIN_0;

        if (context.getScale() > minVisibleScale) {
            Rectangle visibleRect = context.getVisibleRect().intersection(rect);
            Graphics2D g = context.getGraphic2DForColor(Color.gray);
            GraphicUtils.drawCenteredText("Zoom in to see alignments.", visibleRect, g);
            return;
        }

        downsampleRect = new Rectangle(rect);
        downsampleRect.height = DOWNAMPLED_ROW_HEIGHT;
        renderDownsampledIntervals(context, downsampleRect);

        alignmentsRect = new Rectangle(rect);
        alignmentsRect.y += DOWNAMPLED_ROW_HEIGHT + DS_MARGIN_2;
        renderAlignments(context, alignmentsRect);
    }

    private void renderDownsampledIntervals(RenderContext context, Rectangle downsampleRect) {

        // Might be offscreen
        if (!context.getVisibleRect().intersects(downsampleRect)) return;

        final Collection<AlignmentInterval> loadedIntervals = dataManager.getLoadedIntervals(context.getReferenceFrame());
        if (loadedIntervals == null) return;

        Graphics2D g = context.getGraphic2DForColor(Color.black);

        for (AlignmentInterval loadedInterval : loadedIntervals) {
            List<DownsampledInterval> intervals = loadedInterval.getDownsampledIntervals();
            for (DownsampledInterval interval : intervals) {
                int x0 = context.bpToScreenPixel(interval.getStart());
                int x1 = context.bpToScreenPixel(interval.getEnd());
                int w = Math.max(1, x1 - x0);
                // If there is room, leave a gap on one side
                if (w > 5) w--;
                // Greyscale from 0 -> 100 downsampled
                //int gray = 200 - interval.getCount();
                //Color color = (gray <= 0 ? Color.black : ColorUtilities.getGrayscaleColor(gray));
                g.fillRect(x0, downsampleRect.y, w, downsampleRect.height);
            }
        }
    }

    private void renderAlignments(RenderContext context, Rectangle inputRect) {

        //log.debug("Render features");
        Map<String, List<AlignmentInterval.Row>> groups = dataManager.getGroups(context, renderOptions);

        Map<String, PEStats> peStats = dataManager.getPEStats();
        if (peStats != null) {
            renderOptions.peStats = peStats;
        }

        if (groups == null) {
            return;
        }


        Rectangle visibleRect = context.getVisibleRect();
        final boolean leaveMargin = getDisplayMode() == DisplayMode.EXPANDED;

        if (renderOptions.isPairedArcView()) {
            maximumHeight = (int) inputRect.getHeight();
            AlignmentRenderer.getInstance().clearCurveMaps();
        } else {
            maximumHeight = Integer.MAX_VALUE;
        }

        // Divide rectangle into equal height levels
        double y = inputRect.getY();
        double h = expandedHeight;
        if (getDisplayMode() != DisplayMode.EXPANDED) {
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
        for (Map.Entry<String, List<AlignmentInterval.Row>> entry : groups.entrySet()) {
            String group = entry.getKey();
            groupNumber++;

            // Loop through the alignment rows for this group
            List<AlignmentInterval.Row> rows = entry.getValue();
            for (AlignmentInterval.Row row : rows) {

                if ((visibleRect != null && y > visibleRect.getMaxY())) {
                    return;
                }
                if (renderOptions.isPairedArcView()) {
                    y = Math.min(getY() + getHeight(), visibleRect.getMaxY());
                    y -= h;
                }

                if (y + h > visibleRect.getY()) {
                    Rectangle rowRectangle = new Rectangle(inputRect.x, (int) y, inputRect.width, (int) h);
                    renderer.renderAlignments(row.alignments, context, rowRectangle,
                            inputRect, renderOptions, leaveMargin, selectedReadNames);
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
    }

    public void clearCaches() {
        dataManager.clear();
    }

    /**
     * Sort alignment rows based on alignments that intersent location
     */
    public void sortRows(SortOption option, ReferenceFrame referenceFrame, double location, String tag) {
        dataManager.sortRows(option, referenceFrame, location, tag);
    }

    /**
     * Sort alignment rows such that alignments that intersect from the center
     * appear left to right by start position
     */
    public void groupAlignments(GroupOption option, ReferenceFrame referenceFrame) {
        if (renderOptions.groupByOption != option) {
            renderOptions.groupByOption = (option == GroupOption.NONE ? null : option);
            dataManager.repackAlignments(referenceFrame, renderOptions);
        }
    }

    public void packAlignments(ReferenceFrame referenceFrame) {
        dataManager.repackAlignments(referenceFrame, renderOptions);
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
     * Copy the contents of the popup text to the system clipboard.
     */
    public void copyFlowSignalDistribution(final TrackClickEvent e, int location) {
        ArrayList<FlowDistribution> dists = getFlowSignalDistribution(e.getFrame(), location);
        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        String json = "";
        for (FlowDistribution dist : dists) {
            json += dist.toJson() + "\n";
        }
        clipboard.setContents(new StringSelection(json), null);
    }

    /**
     * by default, returns both for forward and backward strand
     */
    private ArrayList<FlowDistribution> getFlowSignalDistribution(ReferenceFrame frame, int location) {
        return getFlowSignalDistribution(frame, location, true, true);
    }

    private ArrayList<FlowDistribution> getFlowSignalDistribution(ReferenceFrame frame, int location, boolean forward, boolean reverse) {
        // one for each base!
        ArrayList<TreeMap<Short, Integer>> alleletrees = new ArrayList<TreeMap<Short, Integer>>();

        int nrflows = 0;
        ArrayList<FlowDistribution> alleledist = new ArrayList<FlowDistribution>();
        String bases = "";
        // also store information on read and position


        ArrayList<ArrayList<ReadInfo>> allelereadinfos = new ArrayList<ArrayList<ReadInfo>>();
        for (AlignmentInterval interval : dataManager.getLoadedIntervals()) {
            Iterator<Alignment> alignmentIterator = interval.getAlignmentIterator();
            while (alignmentIterator.hasNext()) {
                Alignment alignment = alignmentIterator.next();
                if ((alignment.isNegativeStrand() && !reverse) || (!alignment.isNegativeStrand() && !forward)) {
                    continue;
                }
                if (!alignment.contains(location)) {
                    continue;
                }
                // we don't want the beginning or the end of the alignment! HP might might give misleading results
                if (alignment.getAlignmentStart() == location || alignment.getAlignmentEnd() == location) {
                    log.info(location + " for read " + alignment.getReadName() + " is at an end, not taking it");
                    continue;
                }
                // also throw away positions near the end if we have the same base until the end if the user preference is set that way
                boolean hideFirstHPs = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.IONTORRENT_FLOWDIST_HIDE_FIRST_HP);

                if (hideFirstHPs) {
                    char baseatpos = (char) alignment.getBase(location);
                    boolean hp = true;
                    for (int pos = alignment.getAlignmentStart(); pos < location; pos++) {
                        if ((char) alignment.getBase(pos) != baseatpos) {
                            hp = false;
                            break;
                        }
                    }
                    if (hp) {
                        log.info("Got all same bases " + baseatpos + " for read " + alignment.getReadName() + " at START.");
                        continue;
                    }
                    hp = true;
                    for (int pos = location + 1; pos < alignment.getAlignmentEnd(); pos++) {
                        if ((char) alignment.getBase(pos) != baseatpos) {
                            hp = false;
                            break;
                        }
                    }
                    if (hp) {
                        log.info("Got all same bases " + baseatpos + " for read " + alignment.getReadName() + " at END");
                        continue;
                    }
                }
                AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
                for (int i = 0; i < blocks.length; i++) {
                    AlignmentBlock block = blocks[i];
                    int posinblock = (int) location - block.getStart();
                    if (!block.contains((int) location) || !block.hasFlowSignals()) {
                        continue;
                    }

                    int flownr = block.getFlowSignalSubContext(posinblock).getFlowOrderIndex();
                    nrflows++;
                    short flowSignal = block.getFlowSignalSubContext(posinblock).getCurrentSignal();

                    char base = (char) block.getBase(posinblock);

                    int whichbase = bases.indexOf(base);
                    TreeMap<Short, Integer> map = null;
                    ArrayList<ReadInfo> readinfos = null;
                    if (whichbase < 0) {
                        bases += base;
                        map = new TreeMap<Short, Integer>();
                        alleletrees.add(map);
                        readinfos = new ArrayList<ReadInfo>();
                        allelereadinfos.add(readinfos);
                    } else {
                        map = alleletrees.get(whichbase);
                        readinfos = allelereadinfos.get(whichbase);
                    }
                    ReadInfo readinfo = new ReadInfo(alignment.getReadName(), flownr, flowSignal, base);
                    readinfos.add(readinfo);
                    if (map.containsKey(flowSignal)) {
                        // increment
                        map.put(flowSignal, map.get(flowSignal) + 1);
                    } else {
                        // insert
                        map.put(flowSignal, 1);
                    }
                }
            }
        }

        String locus = Locus.getFormattedLocusString(frame.getChrName(), (int) location, (int) location);

        int which = 0;
        for (TreeMap<Short, Integer> map : alleletrees) {
            String name = "";
            if (forward && reverse) {
                name += "both strand";
            } else if (forward) {
                name += "forward strand";
            } else {
                name += "reverse strand";
            }
            char base = bases.charAt(which);
            name += ", " + base + ", " + nrflows + " flows";
            String info = locus + ", " + bases;

            FlowDistribution dist = new FlowDistribution(location, nrflows, map, name, base, forward, reverse, info);
            dist.setReadInfos(allelereadinfos.get(which));
            alleledist.add(dist);
            which++;
        }
        return alleledist;
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
     * Split the screen so the current view and mate region are side by side.
     * Need a better name for this method.
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

    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0.0f;
    }

    AlignmentDataManager getDataManager() {
        return dataManager;
    }

    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        if (downsampleRect != null && y > downsampleRect.y && y <= downsampleRect.y + downsampleRect.height) {
            Collection<AlignmentInterval> loadedIntervals = dataManager.getLoadedIntervals(frame);
            if (loadedIntervals == null) {
                return null;
            } else {
                for (AlignmentInterval loadedInterval : loadedIntervals) {
                    List<DownsampledInterval> intervals = loadedInterval.getDownsampledIntervals();
                    DownsampledInterval interval = (DownsampledInterval) FeatureUtils.getFeatureAt(position, 0, intervals);
                    if (interval != null) return interval.getValueString();
                }
                return null;
            }
        } else if (renderOptions.isPairedArcView()) {
            Alignment feature = null;
            //All alignments stacked at the bottom
            double xloc = (position - frame.getOrigin()) / frame.getScale();
            SortedSet<Shape> arcs = AlignmentRenderer.getInstance().curveOverlap(xloc);
            int halfLength = 2;
            int sideLength = 2 * halfLength;
            for (Shape curve : arcs) {
                if (curve.intersects(xloc - halfLength, y - halfLength, sideLength, sideLength)) {
                    feature = AlignmentRenderer.getInstance().getAlignmentForCurve(curve);
                    break;
                }
            }
            return feature == null ? null : feature.getValueString(position, getWindowFunction());
        } else {
            Alignment feature = getAlignmentAt(position, y, frame);
            return feature == null ? null : feature.getValueString(position, getWindowFunction());

        }
    }

    private Alignment getAlignmentAt(double position, int y, ReferenceFrame frame) {

        if (alignmentsRect == null) {
            return null;   // <= not loaded yet
        }
        Map<String, List<AlignmentInterval.Row>> groups = dataManager.getGroupedAlignmentsContaining(position, frame);

        if (groups == null || groups.isEmpty()) {
            return null;
        }

        int h = getDisplayMode() == DisplayMode.EXPANDED ? expandedHeight : squishedHeight;
//        int levelNumber = (y - renderedRect.y) / h;
//        if (levelNumber < 0) {
//            return null;
//        }

        int startY = alignmentsRect.y;
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
                final boolean showJunctions = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_JUNCTION_TRACK);
                dataManager.setShowSpliceJunctions(showJunctions);
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
                if (parent != null) {
                    parent.repaint();
                }
                return true;
            }

        }
        if (IGV.getInstance().isShowDetailsOnClick()) {
            openTooltipWindow(te);
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
        Color c = alignment.isPaired() && alignment.getMate() != null && alignment.getMate().isMapped()
                ? ColorUtilities.randomColor(selectionColorIndex++) : Color.black;
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

        if (dataManager.isShowSpliceJunctions()) {
            attrs.put("showSpliceJunctions", String.valueOf(dataManager.isShowSpliceJunctions()));
        }

        return attrs;
    }

    @Override
    public void restorePersistentState(Map<String, String> attributes) {
        super.restorePersistentState(attributes);    //To change body of overridden methods use File | Settings | File Templates.
        renderOptions.restorePersistentState(attributes);

        String spliceJunctionOption = attributes.get("showSpliceJunctions");
        if (spliceJunctionOption != null) {
            try {
                dataManager.setShowSpliceJunctions(Boolean.parseBoolean(spliceJunctionOption));
            } catch (Exception e) {
                log.error("Error restoring splice junction option: " + spliceJunctionOption);
            }
        }

    }

    public void setViewAsPairs(boolean vAP) {
        // TODO -- generalize this test to all incompatible pairings
        if (vAP && renderOptions.groupByOption == GroupOption.STRAND) {
            boolean ungroup = MessageUtils.confirm("\"View as pairs\" is incompatible with \"Group by strand\". Ungroup?");
            if (ungroup) {
                renderOptions.groupByOption = null;
            } else {
                return;
            }
        }

        dataManager.setViewAsPairs(vAP, renderOptions);
        refresh();
    }

    public boolean isPairedArcView() {
        return this.renderOptions.isPairedArcView();
    }

    public void setPairedArcView(boolean option) {
        if (option == this.isPairedArcView()) {
            return;
        }

        //TODO This is dumb and bad UI design
        //Should use a combo box or something
        if (option) {
            setViewAsPairs(false);
        }

        renderOptions.setPairedArcView(option);
        for (ReferenceFrame frame : FrameManager.getFrames()) {
            dataManager.repackAlignments(frame, renderOptions);
        }
        refresh();
    }

    public static class RenderOptions {

        ShadeBasesOption shadeBasesOption;
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
        private boolean pairedArcView = false;
        public boolean flagZeroQualityAlignments = true;
        Map<String, PEStats> peStats;
        private String colorByTag;
        private String groupByTag;
        private String sortByTag;

        RenderOptions() {
            PreferenceManager prefs = PreferenceManager.getInstance();

            String shadeOptionString = prefs.get(PreferenceManager.SAM_SHADE_BASES);
            if (shadeOptionString.equals("false")) {
                shadeBasesOption = ShadeBasesOption.NONE;
            } else if (shadeOptionString.equals("true")) {
                shadeBasesOption = ShadeBasesOption.QUALITY;
            } else {
                shadeBasesOption = ShadeBasesOption.valueOf(shadeOptionString);
            }
            shadeCenters = prefs.getAsBoolean(PreferenceManager.SAM_SHADE_CENTER);
            flagUnmappedPairs = prefs.getAsBoolean(PreferenceManager.SAM_FLAG_UNMAPPED_PAIR);
            computeIsizes = prefs.getAsBoolean(PreferenceManager.SAM_COMPUTE_ISIZES);
            minInsertSize = prefs.getAsInt(PreferenceManager.SAM_MIN_INSERT_SIZE_THRESHOLD);
            maxInsertSize = prefs.getAsInt(PreferenceManager.SAM_MAX_INSERT_SIZE_THRESHOLD);
            minInsertSizePercentile = prefs.getAsFloat(PreferenceManager.SAM_MIN_INSERT_SIZE_PERCENTILE);
            maxInsertSizePercentile = prefs.getAsFloat(PreferenceManager.SAM_MAX_INSERT_SIZE_PERCENTILE);
            showAllBases = DEFAULT_SHOWALLBASES;
            colorOption = CollUtils.valueOf(ColorOption.class, prefs.get(PreferenceManager.SAM_COLOR_BY), ColorOption.NONE);
            groupByOption = null;
            flagZeroQualityAlignments = prefs.getAsBoolean(PreferenceManager.SAM_FLAG_ZERO_QUALITY);
            bisulfiteContext = DEFAULT_BISULFITE_CONTEXT;


            colorByTag = prefs.get(PreferenceManager.SAM_COLOR_BY_TAG);
            sortByTag = prefs.get(PreferenceManager.SAM_SORT_BY_TAG);
            groupByTag = prefs.get(PreferenceManager.SAM_GROUP_BY_TAG);

            //updateColorScale();

            peStats = new HashMap<String, PEStats>();
        }

        /*
         * private void updateColorScale() { int delta = 1; if (medianInsertSize
         * == 0 || madInsertSize == 0) { delta = (maxInsertSizeThreshold -
         * minInsertSizeThreshold) / 10; } else { delta =
         * Math.min((maxInsertSizeThreshold - minInsertSizeThreshold) / 3,
         * madInsertSize); } insertSizeColorScale = new
         * ContinuousColorScale(minInsertSizeThreshold, minInsertSizeThreshold +
         * delta, maxInsertSizeThreshold - delta, maxInsertSizeThreshold,
         * Color.blue, AlignmentRenderer.grey1, Color.red); }
         */

        /**
         * Called by session writer. Return instance variable values as a map of
         * strings. Used to record current state of object. Variables with
         * default values are not stored, as it is presumed the user has not
         * changed them.
         *
         * @return
         */
        public Map<String, String> getPersistentState() {
            Map<String, String> attributes = new HashMap();
            PreferenceManager prefs = PreferenceManager.getInstance();
            if (shadeBasesOption != DEFAULT_SHADE_BASES_OPTION) {
                attributes.put("shadeBasesOption", shadeBasesOption.toString());
            }
            if (shadeCenters != prefs.getAsBoolean(PreferenceManager.SAM_SHADE_CENTER)) {
                attributes.put("shadeCenters", String.valueOf(shadeCenters));
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
            value = attributes.get("shadeBases");  // For older sessions
            if (value != null) {
                if (value.equals("false")) {
                    shadeBasesOption = ShadeBasesOption.NONE;
                } else if (value.equals("true")) {
                    shadeBasesOption = ShadeBasesOption.QUALITY;
                }
            }
            shadeBasesOption = getFromMap(attributes, "shadeBasesOption", ShadeBasesOption.class, ShadeBasesOption.QUALITY);
            colorOption = getFromMap(attributes, "colorOption", ColorOption.class, ColorOption.NONE);
            groupByOption = getFromMap(attributes, "groupByOption", GroupOption.class, GroupOption.NONE);
            bisulfiteContext = getFromMap(attributes, "bisulfiteContextRenderOption", BisulfiteContext.class, BisulfiteContext.CG);


            groupByTag = getFromMap(attributes, "groupByTag", groupByTag);
            sortByTag = getFromMap(attributes, "sortByTag", sortByTag);
            colorByTag = getFromMap(attributes, "colorByTag", colorByTag);
        }

        private <T extends Enum<T>> T getFromMap(Map<String, String> attributes, String key, Class<T> clazz, T defaultValue) {
            String value = attributes.get(key);
            if (value == null) {
                return defaultValue;
            }
            return CollUtils.<T>valueOf(clazz, value, defaultValue);
        }

        private String getFromMap(Map<String, String> attributes, String key, String defaultValue) {
            String value = attributes.get(key);
            if (value == null) {
                return defaultValue;
            }
            return value;
        }

        public boolean isPairedArcView() {
            return pairedArcView;
        }

        public void setPairedArcView(boolean pairedArcView) {
            this.pairedArcView = pairedArcView;
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

            if (ionTorrent) {
                addCopyFlowSignalDistributionToClipboardItem(e);
                addIonTorrentAuxiliaryViews(e);
            }

            addSeparator();
            addGroupMenuItem();
            addSortMenuItem();
            addColorByMenuItem();

            addSeparator();
            addShadeBaseByMenuItem();
            addShowMismatchesMenuItem();
            addShowAllBasesMenuItem();

            addSeparator();
            addViewAsPairsMenuItem();

            boolean viewPairArcsPresent = Boolean.parseBoolean(System.getProperty("pairedArcViewPresent", "false"));
            if (viewPairArcsPresent) {
                addViewPairedArcsMenuItem();
            }

            addGoToMate(e);
            showMateRegion(e);
            addInsertSizeMenuItem();

            addSeparator();
            addPackMenuItem();
            //addCoverageDepthMenuItem();
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

        private JCheckBoxMenuItem getGroupMenuItem(String label, final GroupOption option) {
            JCheckBoxMenuItem mi = new JCheckBoxMenuItem(label);
            mi.setSelected(renderOptions.groupByOption == option);
            if (option == GroupOption.NONE) {
                mi.setSelected(renderOptions.groupByOption == null);
            }
            mi.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    IGV.getInstance().groupAlignmentTracks(option);
                    refresh();

                }
            });

            return mi;
        }

        public void addGroupMenuItem() {//ReferenceFrame frame) {
            // Change track height by attribute
            JMenu groupMenu = new JMenu("Group alignments by");
            ButtonGroup group = new ButtonGroup();

            Map<String, GroupOption> mappings = new LinkedHashMap<String, GroupOption>();
            mappings.put("none", GroupOption.NONE);
            mappings.put("read strand", GroupOption.STRAND);
            mappings.put("first-in-pair strand", GroupOption.FIRST_OF_PAIR_STRAND);
            mappings.put("sample", GroupOption.SAMPLE);
            mappings.put("read group", GroupOption.READ_GROUP);
            mappings.put("chromosome of mate", GroupOption.MATE_CHROMOSOME);

            String addExtraStr = System.getProperty("enable.groupby.extras", "false");
            boolean addExtras = Boolean.parseBoolean(addExtraStr);
            if (addExtras) {
                mappings.put("pair orientation", GroupOption.PAIR_ORIENTATION);
            }

            for (Map.Entry<String, GroupOption> el : mappings.entrySet()) {
                JCheckBoxMenuItem mi = getGroupMenuItem(el.getKey(), el.getValue());
                groupMenu.add(mi);
                group.add(mi);
            }

            JCheckBoxMenuItem tagOption = new JCheckBoxMenuItem("tag");
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

            add(groupMenu);
        }

        private JMenuItem getSortMenuItem(String label, final SortOption option) {
            JMenuItem mi = new JMenuItem(label);
            mi.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    IGV.getInstance().sortAlignmentTracks(option, null);
                    refresh();

                }
            });

            return mi;
        }

        /**
         * Sort menu
         */
        public void addSortMenuItem() {


            JMenu sortMenu = new JMenu("Sort alignments by");
            //LinkedHashMap is supposed to preserve order of insertion for iteration
            Map<String, SortOption> mappings = new LinkedHashMap<String, SortOption>();

            mappings.put("start location", SortOption.START);
            mappings.put("read strand", SortOption.STRAND);
            mappings.put("first-of-pair strand", SortOption.FIRST_OF_PAIR_STRAND);
            mappings.put("base", SortOption.NUCELOTIDE);
            mappings.put("mapping quality", SortOption.QUALITY);
            mappings.put("sample", SortOption.SAMPLE);
            mappings.put("read group", SortOption.READ_GROUP);

            if (dataManager.isPairedEnd()) {
                mappings.put("insert size", SortOption.INSERT_SIZE);
                mappings.put("chromosome of mate", SortOption.MATE_CHR);
            }

            for (Map.Entry<String, SortOption> el : mappings.entrySet()) {
                sortMenu.add(getSortMenuItem(el.getKey(), el.getValue()));
            }

            JMenuItem tagOption = new JMenuItem("tag");
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

            // TODO Setting "color-by bisulfite"  also controls the experiment type.  This is temporary, until we
            // expose experimentType directory.
            ExperimentType t = (option == ColorOption.BISULFITE ? ExperimentType.BISULFITE : ExperimentType.OTHER);
            setExperimentType(t);

        }

        private JRadioButtonMenuItem getColorMenuItem(String label, final ColorOption option) {
            JRadioButtonMenuItem mi = new JRadioButtonMenuItem(label);
            mi.setSelected(renderOptions.colorOption == option);
            mi.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    setColorOption(option);
                    refresh();
                }
            });

            return mi;
        }

        public void addColorByMenuItem() {
            // Change track height by attribute
            JMenu colorMenu = new JMenu("Color alignments by");

            ButtonGroup group = new ButtonGroup();

            Map<String, ColorOption> mappings = new LinkedHashMap<String, ColorOption>();

            mappings.put("no color", ColorOption.NONE);

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


            for (Map.Entry<String, ColorOption> el : mappings.entrySet()) {
                JRadioButtonMenuItem mi = getColorMenuItem(el.getKey(), el.getValue());
                colorMenu.add(mi);
                group.add(mi);
            }

            JRadioButtonMenuItem tagOption = new JRadioButtonMenuItem("tag");
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
                final Alignment alignment = getAlignmentAt(location, me.getY(), frame);

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

        public void addCopyFlowSignalDistributionToClipboardItem(final TrackClickEvent te) {
            final MouseEvent me = te.getMouseEvent();
            JMenuItem item = new JMenuItem("Copy the flow signal distrubtion for this base to the clipboard");
            final ReferenceFrame frame = te.getFrame();
            if (frame == null) {
                item.setEnabled(false);
            } else {
                final int location = (int) (frame.getChromosomePosition(me.getX()));

                // Change track height by attribute
                item.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {
                        copyFlowSignalDistribution(te, location);
                    }
                });
            }
            add(item);
        }

        public void addViewPairedArcsMenuItem() {
            final JMenuItem item = new JCheckBoxMenuItem("View paired arcs");
            item.setSelected(isPairedArcView());
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    boolean isPairedArcView = item.isSelected();
                    setPairedArcView(isPairedArcView);
                }
            });
            item.setEnabled(dataManager.isPairedEnd());
            add(item);
        }

        public void addViewAsPairsMenuItem() {
            final JMenuItem item = new JCheckBoxMenuItem("View as pairs");
            item.setSelected(renderOptions.isViewPairs());
            item.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    boolean viewAsPairs = item.isSelected();
                    setViewAsPairs(viewAsPairs);
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
                final Alignment alignment = getAlignmentAt(location, e.getY(), frame);
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
                final Alignment alignment = getAlignmentAt(location, e.getY(), frame);

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

        //        public void addCoverageDepthMenuItem() {
//            // Change track height by attribute
//            final JMenuItem item = new JCheckBoxMenuItem("Set maximum coverage depth ...");
//            item.addActionListener(new ActionListener() {
//
//                public void actionPerformed(ActionEvent aEvt) {
//                    int maxLevels = dataManager.getDownsampleOptions().getMaxReadCount();
//                    String val = MessageUtils.showInputDialog("Maximum coverage depth", String.valueOf(maxLevels));
//                    //Cancel button return null
//                    if (val == null) {
//                        return;
//                    }
//                    try {
//                        int newMaxLevels = Integer.parseInt(val);
//                        if (newMaxLevels != maxLevels) {
//                            setMaxDepth(newMaxLevels);
//                        }
//                    } catch (NumberFormatException ex) {
//                        MessageUtils.showMessage("Insert size must be an integer value: " + val);
//                    }
//
//                }
//            });
//            add(item);
//        }
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

        public void addShadeBaseByMenuItem() {

            if (ionTorrent) {
                JMenu groupMenu = new JMenu("Shade bases by...");
                ButtonGroup group = new ButtonGroup();

                Map<String, ShadeBasesOption> mappings = new LinkedHashMap<String, ShadeBasesOption>();
                mappings.put("none", ShadeBasesOption.NONE);
                mappings.put("quality", ShadeBasesOption.QUALITY);
                mappings.put("read flow signal deviation", ShadeBasesOption.FLOW_SIGNAL_DEVIATION_READ);
                mappings.put("reference flow signal deviation", ShadeBasesOption.FLOW_SIGNAL_DEVIATION_REFERENCE);

                for (Map.Entry<String, ShadeBasesOption> el : mappings.entrySet()) {
                    JCheckBoxMenuItem mi = getShadeBasesMenuItem(el.getKey(), el.getValue());
                    groupMenu.add(mi);
                    group.add(mi);
                }

                add(groupMenu);
            } else {
                final JMenuItem item = new JCheckBoxMenuItem("Shade base by quality");
                item.setSelected(renderOptions.shadeBasesOption == ShadeBasesOption.QUALITY);
                item.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {
                        UIUtilities.invokeOnEventThread(new Runnable() {

                            public void run() {
                                if (item.isSelected()) {
                                    renderOptions.shadeBasesOption = ShadeBasesOption.QUALITY;
                                } else {
                                    renderOptions.shadeBasesOption = ShadeBasesOption.NONE;
                                }
                                refresh();
                            }
                        });
                    }
                });

                add(item);
            }

        }

        private JCheckBoxMenuItem getShadeBasesMenuItem(String label, final ShadeBasesOption option) {
            final JCheckBoxMenuItem mi = new JCheckBoxMenuItem(label);
            mi.setSelected(renderOptions.shadeBasesOption == option);

            if (option == ShadeBasesOption.NONE) {
                mi.setSelected(renderOptions.shadeBasesOption == null || renderOptions.shadeBasesOption == option);
            }
            mi.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent aEvt) {
                    UIUtilities.invokeOnEventThread(new Runnable() {

                        public void run() {
                            if (mi.isSelected()) {
                                renderOptions.shadeBasesOption = option;
                            } else {
                                renderOptions.shadeBasesOption = ShadeBasesOption.NONE;
                            }
                            refresh();
                        }
                    });

                }
            });

            return mi;
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
                                    CoverageDataSource ds = new GobyCountArchiveDataSource(file);
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

        private void addIonTorrentAuxiliaryViews(final TrackClickEvent e) {

            JMenu groupMenu = new JMenu("View flow signal distrubtion for this base for");
            ButtonGroup group = new ButtonGroup();


            JCheckBoxMenuItem itemb = new JCheckBoxMenuItem("forward and reverse (2 data series)");
            groupMenu.add(itemb);
            group.add(itemb);
            itemb.setSelected(true);
            itemb.setFont(itemb.getFont().deriveFont(Font.BOLD));

            JCheckBoxMenuItem item = new JCheckBoxMenuItem("both strands combined");
            groupMenu.add(item);
            group.add(item);

            JCheckBoxMenuItem itemf = new JCheckBoxMenuItem("forward strand only");
            groupMenu.add(itemf);
            group.add(itemf);

            JCheckBoxMenuItem itemr = new JCheckBoxMenuItem("reverse strand only");
            groupMenu.add(itemr);
            group.add(itemr);

            final ReferenceFrame frame = e.getFrame();
            if (frame == null) {
                item.setEnabled(false);
                itemf.setEnabled(false);
                itemr.setEnabled(false);
                itemb.setEnabled(false);
            } else {
                final int location = (int) (frame.getChromosomePosition(e.getMouseEvent().getX()));
                // Change track height by attribute
                item.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {
                        showFlowSignalDistribution(location, e.getFrame(), true, true);
                    }
                });
                itemf.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {
                        showFlowSignalDistribution(location, e.getFrame(), true, false);
                    }
                });
                itemr.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {
                        showFlowSignalDistribution(location, e.getFrame(), false, true);
                    }
                });
                itemb.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent aEvt) {
                        showFlowSignalDistribution(location, e.getFrame(), false, false);
                    }
                });
            }
            add(groupMenu);
        }
    }

    /**
     * if neither forward nor reverse, create 2 charts in one
     */
    private void showFlowSignalDistribution(final int location, final ReferenceFrame frame, final boolean forward, final boolean reverse) {
        FlowDistribution[] distributions = getFlowDistributions(forward, reverse, frame, location);

        final FlowSignalDistributionPanel distributionPanel = new FlowSignalDistributionPanel(distributions);
        LocationListener listener = new LocationListener() {

            @Override
            public void locationChanged(int newLocation) {
                log.info("Got new location from panel: " + newLocation + ", (old location was: " + location + ")");
                FlowDistribution[] newdist = getFlowDistributions(forward, reverse, frame, newLocation);
                distributionPanel.setDistributions(newdist);
                //frame.jumpTo(frame.getChrName(), location, location);

                frame.centerOnLocation(newLocation + 1);
                IGV.repaintPanelsHeadlessSafe();
            }
        };
        distributionPanel.setListener(listener);

        // listen to left/right mouse clicks from panel and navigate accordingly
        SimpleDialog dia = new SimpleDialog("Flow Signal Distribution", distributionPanel, 800, 500);
    }

    private FlowDistribution[] getFlowDistributions(boolean forward, boolean reverse, ReferenceFrame frame, int location) {
        FlowDistribution distributions[] = null;
        if (forward || reverse) {
            ArrayList<FlowDistribution> dists = getFlowSignalDistribution(frame, location, forward, reverse);
            distributions = new FlowDistribution[dists.size()];
            for (int i = 0; i < dists.size(); i++) {
                distributions[i] = dists.get(i);
            }
        } else {
            ArrayList<FlowDistribution> distsf = getFlowSignalDistribution(frame, location, true, false);
            ArrayList<FlowDistribution> distsr = getFlowSignalDistribution(frame, location, false, true);
            distributions = new FlowDistribution[distsf.size() + distsr.size()];
            for (int i = 0; i < distsf.size(); i++) {
                distributions[i] = distsf.get(i);
            }
            for (int i = 0; i < distsr.size(); i++) {
                distributions[i + distsf.size()] = distsr.get(i);
            }

        }
        return distributions;
    }

}
