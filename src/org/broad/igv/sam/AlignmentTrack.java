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
package org.broad.igv.sam;

//~--- non-JDK imports --------------------------------------------------------

import com.jidesoft.swing.JidePopupMenu;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.lists.GeneList;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.session.Session;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.panel.*;
import org.broad.igv.ui.util.FileChooserDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ColorUtilities;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.MouseEvent;
import java.io.File;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class AlignmentTrack extends AbstractTrack implements DragListener {


    public enum SortOption {
        START, STRAND, NUCELOTIDE, QUALITY, SAMPLE, READ_GROUP
    }

    public enum ColorOption {
        INSERT_SIZE, READ_STRAND, FRAGMENT_STRAND, PAIR_ORIENTATION, READ_GROUP;
    }


    public static final int MIN_ALIGNMENT_SPACING = 10;
    static final ColorOption DEFAULT_COLOR_OPTION = ColorOption.INSERT_SIZE;
    static final boolean DEFAULT_SHOWALLBASES = false;

    private static ColorOption colorByOption = null;

    private SequenceTrack sequenceTrack;


    private CoverageTrack coverageTrack;

    private RenderOptions renderOptions;

    private static Logger log = Logger.getLogger(AlignmentTrack.class);
    private int expandedHeight = 14;
    private int collapsedHeight = 2;
    private FeatureRenderer renderer;
    private double minVisibleScale = 25;
    private Rectangle renderedRect;
    private HashMap<String, Color> selectedReadNames = new HashMap();
    private int selectionColorIndex = 0;
    private int minHeight = 100;
    private AlignmentDataManager dataManager;

    public CoverageTrack getCoverageTrack() {
        return coverageTrack;
    }


    public static class RenderOptions {
        boolean shadeBases;
        boolean shadeCenters;
        boolean flagUnmappedPairs;
        boolean showAllBases;
        int minInsertSizeThreshold;
        int insertSizeThreshold;
        ColorOption colorOption;

        RenderOptions() {
            PreferenceManager prefs = PreferenceManager.getInstance();
            shadeBases = prefs.getAsBoolean(PreferenceManager.SAM_SHADE_BASE_QUALITY);
            shadeCenters = prefs.getAsBoolean(PreferenceManager.SAM_SHADE_CENTER);
            flagUnmappedPairs = prefs.getAsBoolean(PreferenceManager.SAM_FLAG_UNMAPPED_PAIR);
            minInsertSizeThreshold = prefs.getAsInt(PreferenceManager.SAM_MIN_INSERT_SIZE_THRESHOLD);
            insertSizeThreshold = prefs.getAsInt(PreferenceManager.SAM_INSERT_SIZE_THRESHOLD);
            showAllBases = DEFAULT_SHOWALLBASES;
            colorOption = colorByOption;
        }

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
            if (insertSizeThreshold != prefs.getAsInt(PreferenceManager.SAM_INSERT_SIZE_THRESHOLD)) {
                attributes.put("insertSizeThreshold", String.valueOf(insertSizeThreshold));
            }
            if (minInsertSizeThreshold != prefs.getAsInt(PreferenceManager.SAM_MIN_INSERT_SIZE_THRESHOLD)) {
                attributes.put("minInsertSizeThreshold", String.valueOf(insertSizeThreshold));
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
                insertSizeThreshold = Integer.parseInt(value);
            }
            value = attributes.get("minInsertSizeThreshold");
            if (value != null) {
                minInsertSizeThreshold = Integer.parseInt(value);
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
        }

    }


    public AlignmentTrack(ResourceLocator locator, AlignmentDataManager dataManager) {
        super(locator);

        //AlignmentQueryReader reader = SamQueryReaderFactory.getReader(locator);
        this.dataManager = dataManager;

        PreferenceManager prefs = PreferenceManager.getInstance();

        float maxRange = prefs.getAsFloat(PreferenceManager.SAM_MAX_VISIBLE_RANGE);
        minVisibleScale = (maxRange * 1000) / 700;

        renderer = new AlignmentRenderer();
        this.setExpanded(true);

        if (prefs.getAsBoolean(PreferenceManager.SAM_SHOW_REF_SEQ)) {
            sequenceTrack = new SequenceTrack("Reference");
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
                }
                catch (Exception e) {
                    log.error("Error setting color option", e);
                    colorByOption = DEFAULT_COLOR_OPTION;

                }
            }
        }
    }

    public void setCoverageTrack(CoverageTrack coverageTrack) {
        this.coverageTrack = coverageTrack;
    }

    public void setRenderer(FeatureRenderer renderer) {
        this.renderer = renderer;
    }

    @Override
    public void setHeight(int preferredHeight) {
        super.setHeight(preferredHeight);
        minHeight = preferredHeight;
    }

    /**
     * Method description
     *
     * @return
     */
    @Override
    public int getHeight() {
        int h = Math.max(minHeight, getNLevels() * (getRowHeight()) + 20);
        return h;
    }

    private int getRowHeight() {
        return isExpanded() ? expandedHeight : collapsedHeight;
    }

    private int getNLevels() {
        return dataManager.getNLevels();
    }

    @Override
    public int getPreferredHeight() {
        return Math.max(100, getHeight());
    }


    public void render(RenderContext context, Rectangle rect) {


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

            if (tmp == null) {
                return;
            }


            Rectangle visibleRect = context.getVisibleRect();

            // Divide rectangle into equal height levels
            double y = inputRect.getY();
            double h = expandedHeight;
            if (!isExpanded()) {
                int visHeight = context.getVisibleRect().height;
                collapsedHeight = Math.max(1, Math.min(expandedHeight, visHeight / dataManager.getMaxDepth(context.getReferenceFrame())));
                h = collapsedHeight;
            }

            int levelNumber = 0;
            for (AlignmentInterval.Row row : tmp) {

                if ((visibleRect != null && y > visibleRect.getMaxY()) || levelNumber > dataManager.getMaxLevels()) {
                    return;
                }

                if (y + h > visibleRect.getY()) {
                    Rectangle rect = new Rectangle(inputRect.x, (int) y, inputRect.width, (int) h);
                    renderOptions.colorOption = colorByOption;
                    renderer.renderAlignments(row.alignments,
                            context,
                            rect,
                            renderOptions,
                            isExpanded(),
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
     * Sort alignment rows such that alignments that intersect from the
     * center appear left to right by start position
     */
    public void sortRows(SortOption option, ReferenceFrame referenceFrame) {
        dataManager.sortRows(option, referenceFrame);
    }

    public void sortRows(SortOption option, ReferenceFrame referenceFrame, double location) {
        dataManager.sortRows(option, referenceFrame, location);
    }


    public void packAlignments(ReferenceFrame referenceFrame) {
        dataManager.repackAlignments(referenceFrame);
    }

    /**
     * Copy the contents of the popup text to the system clipboard.
     */
    public void copyToClipboard(final TrackClickEvent e) {
        double location = e.getFrame().getChromosomePosition(e.getMouseEvent().getX());
        double displayLocation = location + 1;
        Alignment alignment = this.getAlignmentAt(displayLocation, e.getMouseEvent().getY(), e.getFrame());

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
    public void gotoMate(final TrackClickEvent te) {

        MouseEvent e = te.getMouseEvent();
        double location = te.getFrame().getChromosomePosition(e.getX());
        double displayLocation = location + 1;
        Alignment alignment = this.getAlignmentAt(displayLocation, e.getY(), te.getFrame());

        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {
                String chr = mate.mateChr;
                int start = mate.mateStart - 1;
                te.getFrame().centerOnLocation(chr, start);
                te.getFrame().recordHistory();
            }
        }
    }

    /**
     * Split the screen so the current view and mate region are side by side.  Need a better
     * name for this method.
     */
    public void splitScreenMate(final TrackClickEvent te) {

        MouseEvent e = te.getMouseEvent();
        double location = te.getFrame().getChromosomePosition(e.getX());
        double displayLocation = location + 1;
        Alignment alignment = this.getAlignmentAt(displayLocation, e.getY(), te.getFrame());

        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {
                String mateChr = mate.mateChr;
                int mateStart = mate.mateStart - 1;

                ReferenceFrame frame = te.getFrame();
                String locus1 = frame.getCurrentLocusString();

                // Generate a locus string for the read mate.  Keep the window width (in base pairs) == to the current range
                ReferenceFrame.Range range = frame.getCurrentRange();
                int length = range.getLength();
                int s2 = Math.max(0, mateStart - length / 2);
                int e2 = s2 + length;
                String startStr = NumberFormat.getInstance().format(s2);
                String endStr = NumberFormat.getInstance().format(e2);
                String mateLocus = mateChr + ":" + startStr + "-" + endStr;

                Session currentSession = IGVMainFrame.getInstance().getSession();

                // If we are already in gene list mode add the mate as another panel, otherwise switch to gl mode
                if (FrameManager.isGeneListMode()) {
                    currentSession.addGene(mateLocus);
                } else {
                    String listName = locus1 + " <-> " + mateLocus;
                    GeneList geneList = new GeneList(listName, Arrays.asList(locus1, mateLocus));
                    currentSession.setCurrentGeneList(geneList);
                }
                IGVMainFrame.getInstance().resetFrames();

            }
        }
    }


    public void setStatType(WindowFunction type) {
        // ignored
    }

    public WindowFunction getWindowFunction() {
        return null;
    }

    public void setRendererClass(Class rc) {
        // ignored
    }

    // SamTracks use posA custom renderer, not derived from Renderer

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

        // TODO -- highlight mate

        String tmp = (feature == null) ? null : feature.getValueString(position, getWindowFunction());

        return tmp;
    }

    private Alignment getAlignmentAt(double position, int y, ReferenceFrame frame) {

        List<AlignmentInterval.Row> alignmentRows = dataManager.getAlignmentRows(frame);

        if (alignmentRows == null || alignmentRows.isEmpty()) {
            return null;
        }

        int h = isExpanded() ? expandedHeight : collapsedHeight;
        int levelNumber = (y - renderedRect.y) / h;
        if (levelNumber < 0 || levelNumber >= alignmentRows.size()) {
            return null;
        }

        AlignmentInterval.Row row = alignmentRows.get(levelNumber);
        List<Alignment> features = row.alignments;

        // give posA 2 pixel window, otherwise very narrow features will be missed.
        double bpPerPixel = frame.getScale();
        double minWidth = 2 * bpPerPixel;    /* * */
        return (Alignment) FeatureUtils.getFeatureAt(position, minWidth, features);

    }

    public void dragStopped(DragEvent evt) {
        // Disabled.  Not sure why we ever thought this was posA good idea
        //if (PreferenceManager.getInstance().getSAMPreferences().isAutosort() &&
        //        ReferenceFrame.getInstance().getScale() < 1) {
        //    sortRows(SortOption.START);
        //    IGVMainFrame.getInstance().repaintDataPanels();
        //}
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
    public boolean handleClick(TrackClickEvent te) {
        MouseEvent e = te.getMouseEvent();
        if (e.isPopupTrigger()) {
            getPopupMenu(te).show(e.getComponent(), e.getX(), e.getY());
            return true;
        } else {
            if (e.isShiftDown() || e.isAltDown() || (e.getClickCount() > 1)) {
                return super.handleClick(te);
            } else if (e.getButton() == MouseEvent.BUTTON1 &&
                    (Globals.IS_MAC && e.isMetaDown() || (!Globals.IS_MAC && e.isControlDown()))) {
                double location = te.getFrame().getChromosomePosition(e.getX());
                double displayLocation = location + 1;
                Alignment alignment = this.getAlignmentAt(displayLocation, e.getY(), te.getFrame());
                if (alignment != null) {
                    if (selectedReadNames.containsKey(alignment.getReadName())) {
                        selectedReadNames.remove(alignment.getReadName());
                    } else {
                        Color c = alignment.isPaired() && alignment.getMate() != null && alignment.getMate().isMapped() ?
                                ColorUtilities.randomColor(selectionColorIndex++) : Color.black;
                        selectedReadNames.put(alignment.getReadName(), c);
                    }
                    Object source = e.getSource();
                    if (source instanceof JComponent) {
                        ((JComponent) source).repaint();
                    }
                }
                return true;
            }
        }
        return false;
    }

    public JPopupMenu getPopupMenu(final TrackClickEvent e) {

        MouseEvent evt = e.getMouseEvent();

        JPopupMenu popupMenu = new JidePopupMenu();

        JLabel popupTitle = new JLabel("  " + getName(), JLabel.CENTER);

        Font newFont = popupMenu.getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        if (popupTitle != null) {
            popupMenu.add(popupTitle);
        }

        addSortMenuItem(popupMenu, e.getFrame());
        addPackMenuItem(popupMenu, e.getFrame());
        addCoverageDepthMenuItem(popupMenu);
        popupMenu.addSeparator();
        addColorByMenuItem(popupMenu);

        addShadeBaseMenuItem(popupMenu);
        addShadeCentersMenuItem(popupMenu);

        addGoToMate(popupMenu, e);
        showMateRegion(popupMenu, e);
        addShowAllBasesMenuItem(popupMenu);
        addMinInsertSizeMenuItem(popupMenu);
        addMaxInsertSizeMenuItem(popupMenu);
        popupMenu.addSeparator();

        addShowCoverageItem(popupMenu);
        addLoadCoverageDataItem(popupMenu);

        addCopyToClipboardItem(popupMenu, e);
        popupMenu.addSeparator();

        addSelecteByNameItem(popupMenu);
        popupMenu.addSeparator();

        JLabel trackSettingsHeading = new JLabel("  Track Settings", JLabel.LEFT);
        trackSettingsHeading.setFont(newFont);

        popupMenu.add(trackSettingsHeading);


        Collection<Track> tmp = new ArrayList();
        tmp.add(this);
        popupMenu.add(TrackMenuUtils.getTrackRenameItem(tmp));

        popupMenu.add(TrackMenuUtils.getExpandCollapseItem(tmp));

        popupMenu.add(TrackMenuUtils.getRemoveMenuItem(tmp));

        popupMenu.addSeparator();
        addClearSelectionsMenuItem(popupMenu);

        return popupMenu;
    }

    public void addSelecteByNameItem(JPopupMenu menu) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Select by name...");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                String val = MessageUtils.showInputDialog("Enter read name: ");
                if (val != null && val.trim().length() > 0) {
                    selectedReadNames.put(val, ColorUtilities.randomColor(selectedReadNames.size() + 1));
                    refresh();
                }
            }
        });

        menu.add(item);
    }

    public void addSortMenuItem(JPopupMenu menu, final ReferenceFrame frame) {
        // Change track height by attribute
        JMenu item = new JMenu("Sort alignments");

        JMenuItem m1 = new JMenuItem("by start location");
        m1.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                IGVMainFrame.getInstance().getTrackManager().sortAlignmentTracks(SortOption.START, frame);
                refresh();

            }
        });

        JMenuItem m2 = new JMenuItem("by strand");
        m2.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                IGVMainFrame.getInstance().getTrackManager().sortAlignmentTracks(SortOption.STRAND, frame);
                refresh();

            }
        });

        JMenuItem m3 = new JMenuItem("by base");
        m3.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {

                IGVMainFrame.getInstance().getTrackManager().sortAlignmentTracks(SortOption.NUCELOTIDE, frame);
                refresh();

            }
        });

        JMenuItem m4 = new JMenuItem("by mapping quality");
        m4.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {

                IGVMainFrame.getInstance().getTrackManager().sortAlignmentTracks(SortOption.QUALITY, frame);
                refresh();

            }
        });

        JMenuItem m5 = new JMenuItem("by sample");
        m5.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {

                IGVMainFrame.getInstance().getTrackManager().sortAlignmentTracks(SortOption.SAMPLE, frame);
                refresh();

            }
        });

        JMenuItem m6 = new JMenuItem("by read group");
        m6.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {

                IGVMainFrame.getInstance().getTrackManager().sortAlignmentTracks(SortOption.READ_GROUP, frame);
                refresh();

            }
        });

        if (frame != null && frame.getScale() >= MIN_ALIGNMENT_SPACING) {
            item.setEnabled(false);
        }


        item.add(m1);
        item.add(m2);
        item.add(m3);
        item.add(m4);
        item.add(m5);
        item.add(m6);
        menu.add(item);
    }


    private void setColorOption(ColorOption option) {
        colorByOption = option;
        PreferenceManager.getInstance().put(PreferenceManager.SAM_COLOR_BY, option.toString());
    }

    public void addColorByMenuItem(JPopupMenu menu) {
        // Change track height by attribute
        JMenu item = new JMenu("Color alignments");

        ButtonGroup group = new ButtonGroup();

        JRadioButtonMenuItem m1 = new JRadioButtonMenuItem("by insert size");
        m1.setSelected(colorByOption == ColorOption.INSERT_SIZE);
        m1.addActionListener(new TrackMenuUtils.TrackActionListener() {
            public void action() {
                setColorOption(ColorOption.INSERT_SIZE);
                refresh();
            }
        });

        JRadioButtonMenuItem m1a = new JRadioButtonMenuItem("by pair orientation");
        m1a.setSelected(colorByOption == ColorOption.PAIR_ORIENTATION);
        m1a.addActionListener(new TrackMenuUtils.TrackActionListener() {
            public void action() {
                setColorOption(ColorOption.PAIR_ORIENTATION);
                refresh();
            }
        });

        JRadioButtonMenuItem m2 = new JRadioButtonMenuItem("by read strand");
        m2.setSelected(colorByOption == ColorOption.READ_STRAND);
        m2.addActionListener(new TrackMenuUtils.TrackActionListener() {
            public void action() {
                setColorOption(ColorOption.READ_STRAND);
                refresh();
            }
        });

        JRadioButtonMenuItem m3 = new JRadioButtonMenuItem("by first-in-pair read strand");
        m3.setSelected(colorByOption == ColorOption.FRAGMENT_STRAND);
        m3.addActionListener(new TrackMenuUtils.TrackActionListener() {
            public void action() {
                setColorOption(ColorOption.FRAGMENT_STRAND);
                refresh();
            }
        });

        JRadioButtonMenuItem m4 = new JRadioButtonMenuItem("by read group");
        m4.setSelected(colorByOption == ColorOption.READ_GROUP);
        m4.addActionListener(new TrackMenuUtils.TrackActionListener() {
            public void action() {
                setColorOption(ColorOption.READ_GROUP);
                refresh();
            }
        });

        item.add(m1);
        item.add(m1a);
        item.add(m2);
        item.add(m3);
        item.add(m4);
        group.add(m1);
        group.add(m1a);
        group.add(m2);
        group.add(m3);
        group.add(m4);
        menu.add(item);

    }


    public void addPackMenuItem(JPopupMenu menu, final ReferenceFrame frame) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Re-pack alignments");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {
                        IGVMainFrame.getInstance().getTrackManager().packAlignmentTracks(frame);
                        refresh();
                    }
                });
            }
        });

        menu.add(item);
    }

    public void addCopyToClipboardItem(JPopupMenu menu, final TrackClickEvent te) {

        final MouseEvent me = te.getMouseEvent();

        // Change track height by attribute
        JMenuItem item = new JMenuItem("Copy read details to clipboard");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {
                        copyToClipboard(te);
                    }
                });
            }
        });
        if (te.getFrame() == null || te.getFrame().getScale() >= MIN_ALIGNMENT_SPACING) {
            item.setEnabled(false);
        }

        menu.add(item);
    }

    public void addGoToMate(JPopupMenu menu, final TrackClickEvent te) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Go to mate");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {
            public void action() {
                gotoMate(te);
            }
        });
        if (te.getFrame() == null) {
            item.setEnabled(false);
        }
        menu.add(item);
    }


    public void showMateRegion(JPopupMenu menu, final TrackClickEvent te) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Show mate region");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {
            public void action() {
                splitScreenMate(te);
            }
        });
        if (te.getFrame() == null) {
            item.setEnabled(false);
        }
        menu.add(item);
    }

    public void addClearSelectionsMenuItem(JPopupMenu menu) {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Clear selections");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                selectedReadNames.clear();
                refresh();
            }
        });

        menu.add(item);
    }

    public void addShowAllBasesMenuItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Show all bases");
        item.setSelected(renderOptions.showAllBases);
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {
                        renderOptions.showAllBases = item.isSelected();
                        refresh();
                    }
                });
            }
        });

        menu.add(item);
    }

    public void addCoverageDepthMenuItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Set maximum coverage depth ...");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                int maxLevels = dataManager.getMaxLevels();
                String val = MessageUtils.showInputDialog("Maximum coverage depth", String.valueOf(maxLevels));
                try {
                    int newMaxLevels = Integer.parseInt(val);
                    if (newMaxLevels != maxLevels) {
                        dataManager.setMaxLevels(newMaxLevels);
                        //dataManager.reload();
                        refresh();
                    }
                }
                catch (NumberFormatException ex) {
                    MessageUtils.showMessage("Insert size must be an integer value: " + val);
                }

            }
        });


        menu.add(item);
    }

    public void addMaxInsertSizeMenuItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Set maximum insert size threshold ...");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                int threshold = renderOptions.insertSizeThreshold;
                String val = MessageUtils.showInputDialog("maximum insert size threshold", String.valueOf(threshold));
                try {
                    int newThreshold = Integer.parseInt(val);
                    if (newThreshold != threshold) {
                        renderOptions.insertSizeThreshold = newThreshold;
                        refresh();
                    }
                }
                catch (NumberFormatException ex) {
                    MessageUtils.showMessage("Insert size must be an integer value: " + val);
                }

            }
        });


        menu.add(item);
    }

    public void addMinInsertSizeMenuItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Set minimum insert size threshold ...");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                int threshold = renderOptions.minInsertSizeThreshold;
                String val = MessageUtils.showInputDialog("minimum insert size threshold", String.valueOf(threshold));
                try {
                    int newThreshold = Integer.parseInt(val);
                    if (newThreshold != threshold) {
                        renderOptions.minInsertSizeThreshold = newThreshold;
                        refresh();
                    }
                }
                catch (NumberFormatException ex) {
                    MessageUtils.showMessage("Insert size must be an integer value: " + val);
                }

            }
        });


        menu.add(item);
    }

    private void refresh() {
        IGVMainFrame.getInstance().repaintDataPanels();
    }

    public void addShadeBaseMenuItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Shade base by quality");
        item.setSelected(renderOptions.shadeBases);
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {
                        renderOptions.shadeBases = item.isSelected();
                        refresh();
                    }
                });
            }
        });

        menu.add(item);
    }

    public void addShadeCentersMenuItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Shade alignments intersecting center");
        item.setSelected(renderOptions.shadeCenters);
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {

                        renderOptions.shadeCenters = item.isSelected();
                        refresh();
                    }
                });
            }
        });

        menu.add(item);
    }


    public void addShowCoverageItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Show coverage track");
        item.setSelected(getCoverageTrack() != null && getCoverageTrack().isVisible());
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {
                        if (getCoverageTrack() != null) {
                            getCoverageTrack().setVisible(item.isSelected());
                            refresh();
                            IGVMainFrame.getInstance().repaintNamePanels();
                        }
                    }
                });
            }
        });

        menu.add(item);
    }

    public void addLoadCoverageDataItem(JPopupMenu menu) {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Load coverage data...");
        item.addActionListener(new TrackMenuUtils.TrackActionListener() {

            public void action() {
                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {

                        FileChooserDialog trackFileDialog = IGVMainFrame.getInstance().getTrackFileChooser();
                        trackFileDialog.setMultiSelectionEnabled(false);
                        trackFileDialog.setVisible(true);
                        if (!trackFileDialog.isCanceled()) {
                            File file = trackFileDialog.getSelectedFile();
                            String path = file.getAbsolutePath();
                            if (path.endsWith(".tdf") || path.endsWith(".tdf")) {

                                TDFReader reader = TDFReader.getReader(file.getAbsolutePath());
                                TDFDataSource ds = new TDFDataSource(reader, 0, getName() + " coverage");
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

        menu.add(item);
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
}
