package org.broad.igv.sam;

import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Locatable;
import org.broad.igv.Globals;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.Strand;
import org.broad.igv.jbrowse.CircularViewUtilities;
import org.broad.igv.lists.GeneList;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sashimi.SashimiPlot;
import org.broad.igv.tools.PFMExporter;
import org.broad.igv.track.SequenceTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.InsertSizeSettingsDialog;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.util.extview.ExtendViewClient;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import static org.broad.igv.prefs.Constants.*;

/**
 * Popup menu class for AlignmentTrack.  The menu gets instantiated from TrackPanelComponent on right-click in the
 * alignment track or its associated name panel.
 */
class AlignmentTrackMenu extends IGVPopupMenu {

    private static Logger log = LogManager.getLogger(AlignmentTrackMenu.class);

    private final AlignmentTrack alignmentTrack;
    private final AlignmentDataManager dataManager;
    private final AlignmentTrack.RenderOptions renderOptions;
    private static int nClusters = 2;

    AlignmentTrackMenu(AlignmentTrack alignmentTrack, final TrackClickEvent e) {

        this.alignmentTrack = alignmentTrack;
        this.dataManager = alignmentTrack.getDataManager();
        this.renderOptions = alignmentTrack.getRenderOptions();
        final Alignment clickedAlignment = alignmentTrack.getAlignmentAt(e);

        // Title
        JLabel popupTitle = new JLabel("  " + alignmentTrack.getName(), JLabel.CENTER);
        Font newFont = getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        add(popupTitle);

        // Circular view items -- optional
        if (PreferencesManager.getPreferences().getAsBoolean(CIRC_VIEW_ENABLED) && CircularViewUtilities.ping()) {
            addSeparator();
            JMenuItem item = new JMenuItem("Add Discordant Pairs to Circular View");
            item.setEnabled(alignmentTrack.getDataManager().isPairedEnd());
            add(item);
            item.addActionListener(ae -> alignmentTrack.sendPairsToCircularView(e));

            JMenuItem item2 = new JMenuItem("Add Split Reads to Circular View");
            add(item2);
            item2.addActionListener(ae -> alignmentTrack.sendSplitToCircularView(e));
        }

        // Some generic items from TrackMenuUtils
        Collection<Track> tracks = List.of(alignmentTrack);
        addSeparator();
        add(TrackMenuUtils.getTrackRenameItem(tracks));
        addCopyToClipboardItem(e, clickedAlignment);

        addSeparator();
        JMenuItem item = new JMenuItem("Change Track Color...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackColor(tracks));
        add(item);

        // Experiment type  (RNA, THIRD GEN, OTHER)
        addSeparator();
        addExperimentTypeMenuItem();
        if (alignmentTrack.getExperimentType() == AlignmentTrack.ExperimentType.THIRD_GEN) {
            addHaplotype(e);
        }

        // Linked read items
        addLinkedReadItems();

        // Group, sort, color, shade, and pack
        addSeparator();
        addGroupMenuItem(e);
        addSortMenuItem();
        addColorByMenuItem();
        addShadeAlignmentsMenuItem();
        //addFilterMenuItem();
        addPackMenuItem();

        // Shading and mismatch items
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
        addThirdGenItems(clickedAlignment, e);

        // Display mode items
        addSeparator();
        TrackMenuUtils.addDisplayModeItems(tracks, this);

        // Select alignment items
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

        // Insertion items, only if clicked over an insertion
        AlignmentBlock insertion = alignmentTrack.getInsertion(clickedAlignment, e.getMouseEvent().getX());
        if (insertion != null) {
            addSeparator();
            addInsertionItems(insertion);
        }

        // Sashimi plot, probably should be depdenent on experimentType (RNA)
        addSeparator();
        JMenuItem sashimi = new JMenuItem("Sashimi Plot");
        sashimi.addActionListener(e1 -> SashimiPlot.openSashimiPlot());
        add(sashimi);

        // Show alignments, coverage, splice junctions
        addSeparator();
        addShowItems();
    }

    private void addShowChimericRegions(final AlignmentTrack alignmentTrack, final TrackClickEvent e, final Alignment clickedAlignment) {

        JMenuItem item = new JMenuItem("View chimeric alignments in split screen");
        if (clickedAlignment != null && clickedAlignment.getAttribute(SAMTag.SA.name()) != null) {
            item.setEnabled(true);
            item.addActionListener(aEvt -> {
                final String saTag = clickedAlignment.getAttribute(SAMTag.SA.name()).toString();
                try {
                    List<SupplementaryAlignment> supplementaryAlignments = SupplementaryAlignment.parseFromSATag(saTag);
                    alignmentTrack.setSelectedAlignment(clickedAlignment);
                    addNewLociToFrames(e.getFrame(), supplementaryAlignments);
                } catch (final Exception ex) {
                    MessageUtils.showMessage("Failed to handle SA tag: " + saTag + " due to " + ex.getMessage());
                    item.setEnabled(false);
                }
            });
        } else {
            item.setEnabled(false);
        }
        add(item);
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

            String nString = MessageUtils.showInputDialog("Enter the number of clusters", String.valueOf(nClusters));
            if (nString == null) {
                return;
            }
            try {
                nClusters = Integer.parseInt(nString);
            } catch (NumberFormatException e1) {
                MessageUtils.showMessage("Clusters size must be an integer");
                return;
            }

            final int start = (int) frame.getOrigin();
            final int end = (int) frame.getEnd();

            AlignmentInterval interval = dataManager.getLoadedInterval(frame);
            HaplotypeUtils haplotypeUtils = new HaplotypeUtils(interval);
            boolean success = haplotypeUtils.clusterAlignments(frame.getChrName(), start, end, nClusters);

            if (success) {
                groupAlignments(AlignmentTrack.GroupOption.HAPLOTYPE, null, null);
                alignmentTrack.repaint();
            }

            //dataManager.sortRows(SortOption.HAPLOTYPE, frame, (end + start) / 2, null);
            //AlignmentTrack.repaint();

        });


    }


    /**
     * Item for exporting "consensus" sequence of region, based on loaded alignments.
     *
     * @param e
     */
    private void addConsensusSequence(TrackClickEvent e) {

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

            if (frame == null) {  // Should never happen
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

        JMenu bisulfiteContextMenu = new JMenu("bisulfite mode");

        JRadioButtonMenuItem nomeESeqOption = null;
        boolean showNomeESeq = alignmentTrack.getPreferences().getAsBoolean(SAM_NOMESEQ_ENABLED);
        if (showNomeESeq) {
            nomeESeqOption = new JRadioButtonMenuItem("NOMe-seq bisulfite mode");
            nomeESeqOption.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.NOMESEQ);
            nomeESeqOption.addActionListener(aEvt -> {
                alignmentTrack.setColorOption(AlignmentTrack.ColorOption.NOMESEQ);
                alignmentTrack.repaint();
            });
            group.add(nomeESeqOption);
        }

        for (final AlignmentTrack.BisulfiteContext item : AlignmentTrack.BisulfiteContext.values()) {

            JRadioButtonMenuItem m1 = new JRadioButtonMenuItem(item.getLabel());
            m1.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.BISULFITE && renderOptions.bisulfiteContext == item);
            m1.addActionListener(aEvt -> {
                alignmentTrack.setColorOption(AlignmentTrack.ColorOption.BISULFITE);
                alignmentTrack.setBisulfiteContext(item);
                alignmentTrack.repaint();
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
                alignmentTrack.getSelectedReadNames().put(val, alignmentTrack.getReadNamePalette().get(val));
                alignmentTrack.repaint();
            }
        });
        add(item);
    }

    void addExperimentTypeMenuItem() {
        Map<String, AlignmentTrack.ExperimentType> mappings = new LinkedHashMap<>();
        mappings.put("Other", AlignmentTrack.ExperimentType.OTHER);
        mappings.put("RNA", AlignmentTrack.ExperimentType.RNA);
        mappings.put("3rd Gen", AlignmentTrack.ExperimentType.THIRD_GEN);
        //mappings.put("Bisulfite", ExperimentType.BISULFITE);
        JMenu groupMenu = new JMenu("Experiment Type");
        ButtonGroup group = new ButtonGroup();
        for (Map.Entry<String, AlignmentTrack.ExperimentType> el : mappings.entrySet()) {
            JCheckBoxMenuItem mi = getExperimentTypeMenuItem(el.getKey(), el.getValue());
            groupMenu.add(mi);
            group.add(mi);
        }
        add(groupMenu);
    }

    private JCheckBoxMenuItem getExperimentTypeMenuItem(String label, final AlignmentTrack.ExperimentType option) {
        JCheckBoxMenuItem mi = new JCheckBoxMenuItem(label);
        mi.setSelected(alignmentTrack.getExperimentType() == option);
        mi.addActionListener(aEvt -> alignmentTrack.setExperimentType(option));
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
        final int chromStart = (int) frame.getChromosomePosition(me);

        // Change track height by attribute
        JMenu groupMenu = new JMenu("Group alignments by");
        ButtonGroup group = new ButtonGroup();

        AlignmentTrack.GroupOption[] groupOptions = {
                AlignmentTrack.GroupOption.NONE, AlignmentTrack.GroupOption.STRAND, AlignmentTrack.GroupOption.FIRST_OF_PAIR_STRAND, AlignmentTrack.GroupOption.SAMPLE,
                AlignmentTrack.GroupOption.LIBRARY, AlignmentTrack.GroupOption.READ_GROUP, AlignmentTrack.GroupOption.MATE_CHROMOSOME,
                AlignmentTrack.GroupOption.PAIR_ORIENTATION, AlignmentTrack.GroupOption.SUPPLEMENTARY, AlignmentTrack.GroupOption.REFERENCE_CONCORDANCE,
                AlignmentTrack.GroupOption.MOVIE, AlignmentTrack.GroupOption.ZMW, AlignmentTrack.GroupOption.READ_ORDER, AlignmentTrack.GroupOption.LINKED, AlignmentTrack.GroupOption.PHASE,
                AlignmentTrack.GroupOption.MAPPING_QUALITY
        };

        for (final AlignmentTrack.GroupOption option : groupOptions) {
            JCheckBoxMenuItem mi = new JCheckBoxMenuItem(option.label);
            mi.setSelected(renderOptions.getGroupByOption() == option);
            mi.addActionListener(aEvt -> {
                groupAlignments(option, null, null);
            });
            groupMenu.add(mi);
            group.add(mi);
        }

        JCheckBoxMenuItem tagOption = new JCheckBoxMenuItem("tag");
        tagOption.addActionListener(aEvt -> {
            String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getGroupByTag());
            if (tag != null) {
                if (tag.trim().length() > 0) {
                    groupAlignments(AlignmentTrack.GroupOption.TAG, tag, null);
                } else {
                    groupAlignments(AlignmentTrack.GroupOption.NONE, null, null);
                }
            }

        });
        tagOption.setSelected(renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.TAG);
        groupMenu.add(tagOption);
        group.add(tagOption);

        Range oldGroupByPos = renderOptions.getGroupByPos();
        if (oldGroupByPos != null && renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.BASE_AT_POS) { // already sorted by the base at a position
            JCheckBoxMenuItem oldGroupByPosOption = new JCheckBoxMenuItem("base at " + oldGroupByPos.getChr() +
                    ":" + Globals.DECIMAL_FORMAT.format(1 + oldGroupByPos.getStart()));
            groupMenu.add(oldGroupByPosOption);
            oldGroupByPosOption.setSelected(true);
        }

        if (renderOptions.getGroupByOption() != AlignmentTrack.GroupOption.BASE_AT_POS || oldGroupByPos == null ||
                !oldGroupByPos.getChr().equals(chrom) || oldGroupByPos.getStart() != chromStart) { // not already sorted by this position
            JCheckBoxMenuItem newGroupByPosOption = new JCheckBoxMenuItem("base at " + chrom +
                    ":" + Globals.DECIMAL_FORMAT.format(1 + chromStart));
            newGroupByPosOption.addActionListener(aEvt -> {
                Range groupByPos = new Range(chrom, chromStart, chromStart + 1);
                groupAlignments(AlignmentTrack.GroupOption.BASE_AT_POS, null, groupByPos);
            });
            groupMenu.add(newGroupByPosOption);
            group.add(newGroupByPosOption);
        }

        groupMenu.add(new Separator());
        JCheckBoxMenuItem invertGroupNameSortingOption = new JCheckBoxMenuItem("Reverse group order");
        invertGroupNameSortingOption.setSelected(renderOptions.isInvertGroupSorting());
        invertGroupNameSortingOption.addActionListener(aEvt -> {
            renderOptions.setInvertGroupSorting(!renderOptions.isInvertGroupSorting());
            dataManager.packAlignments(renderOptions);
            alignmentTrack.repaint();
        });
        groupMenu.add(invertGroupNameSortingOption);

        JCheckBoxMenuItem groupAllOption = new JCheckBoxMenuItem("Group all tracks");
        groupAllOption.setSelected(alignmentTrack.getPreferences().getAsBoolean(SAM_GROUP_ALL));
        groupAllOption.addActionListener(aEvt -> {
            alignmentTrack.getPreferences().put(SAM_GROUP_ALL, groupAllOption.getState());
        });
        groupMenu.add(groupAllOption);

        add(groupMenu);
    }

    private void groupAlignments(AlignmentTrack.GroupOption option, String tag, Range pos) {

        if(alignmentTrack.getPreferences().getAsBoolean(SAM_GROUP_ALL)) {
            for(AlignmentTrack t : IGV.getInstance().getAlignmentTracks()) {
                t.groupAlignments(option, tag, pos);
            }
        } else {
            alignmentTrack.groupAlignments(option, tag, pos);
        }
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
        mappings.put("aligned read length", SortOption.ALIGNED_READ_LENGTH);
// mappings.put("supplementary flag", SortOption.SUPPLEMENTARY);

        if (dataManager.isPairedEnd()) {
            mappings.put("insert size", SortOption.INSERT_SIZE);
            mappings.put("chromosome of mate", SortOption.MATE_CHR);
        }

        for (Map.Entry<String, SortOption> el : mappings.entrySet()) {
            JMenuItem mi = new JMenuItem(el.getKey());
            mi.addActionListener(aEvt -> {
                final SortOption option = el.getValue();
                renderOptions.setSortOption(option);
                alignmentTrack.sortAlignmentTracks(option, null, renderOptions.isInvertSorting());
            });
            sortMenu.add(mi);
        }

        JMenuItem tagOption = new JMenuItem("tag");
        tagOption.addActionListener(aEvt -> {
            String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getSortByTag());
            if (tag != null && tag.trim().length() > 0) {
                renderOptions.setSortByTag(tag);
                renderOptions.setSortOption((SortOption.TAG));
                alignmentTrack.sortAlignmentTracks(SortOption.TAG, tag, renderOptions.isInvertSorting());
            }
        });
        sortMenu.add(tagOption);

        sortMenu.add(new Separator());
        JCheckBoxMenuItem invertGroupNameSortingOption = new JCheckBoxMenuItem("reverse sorting");
        invertGroupNameSortingOption.setSelected(renderOptions.isInvertSorting());
        invertGroupNameSortingOption.addActionListener(aEvt -> {
            final boolean updatedInvertSorting = !renderOptions.isInvertSorting();
            renderOptions.setInvertSorting(updatedInvertSorting);
            alignmentTrack.sortAlignmentTracks(renderOptions.getSortOption(), renderOptions.getSortByTag(), updatedInvertSorting);
        });
        sortMenu.add(invertGroupNameSortingOption);
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

    private JRadioButtonMenuItem getColorMenuItem(String label, final AlignmentTrack.ColorOption option) {
        JRadioButtonMenuItem mi = new JRadioButtonMenuItem(label);
        mi.setSelected(renderOptions.getColorOption() == option);
        mi.addActionListener(aEvt -> {
            alignmentTrack.setColorOption(option);
            alignmentTrack.repaint();
        });

        return mi;
    }

    void addColorByMenuItem() {
        // Change track height by attribute
        JMenu colorMenu = new JMenu("Color alignments by");

        ButtonGroup group = new ButtonGroup();

        Map<String, AlignmentTrack.ColorOption> mappings = new LinkedHashMap<>();

        mappings.put("none", AlignmentTrack.ColorOption.NONE);

        if (dataManager.hasYCTags()) {
            mappings.put("YC tag", AlignmentTrack.ColorOption.YC_TAG);
        }

        if (dataManager.isPairedEnd()) {
            mappings.put("insert size", AlignmentTrack.ColorOption.INSERT_SIZE);
            mappings.put("pair orientation", AlignmentTrack.ColorOption.PAIR_ORIENTATION);
            mappings.put("insert size and pair orientation", AlignmentTrack.ColorOption.UNEXPECTED_PAIR);
        }

        mappings.put("read strand", AlignmentTrack.ColorOption.READ_STRAND);

        if (dataManager.isPairedEnd()) {
            mappings.put("first-of-pair strand", AlignmentTrack.ColorOption.FIRST_OF_PAIR_STRAND);
        }

        mappings.put("read group", AlignmentTrack.ColorOption.READ_GROUP);

        if (dataManager.isPairedEnd()) {
            mappings.put("read order", AlignmentTrack.ColorOption.READ_ORDER);
        }

        mappings.put("sample", AlignmentTrack.ColorOption.SAMPLE);
        mappings.put("library", AlignmentTrack.ColorOption.LIBRARY);
        mappings.put("movie", AlignmentTrack.ColorOption.MOVIE);
        mappings.put("ZMW", AlignmentTrack.ColorOption.ZMW);

        for (Map.Entry<String, AlignmentTrack.ColorOption> el : mappings.entrySet()) {
            JRadioButtonMenuItem mi = getColorMenuItem(el.getKey(), el.getValue());
            colorMenu.add(mi);
            group.add(mi);
        }

        JRadioButtonMenuItem tagOption = new JRadioButtonMenuItem("tag");
        tagOption.setSelected(renderOptions.getColorOption() == AlignmentTrack.ColorOption.TAG);
        tagOption.addActionListener(aEvt -> {
            alignmentTrack.setColorOption(AlignmentTrack.ColorOption.TAG);
            String tag = MessageUtils.showInputDialog("Enter tag", renderOptions.getColorByTag());
            if (tag != null && tag.trim().length() > 0) {
                alignmentTrack.setColorByTag(tag);
                alignmentTrack.repaint();
            }
        });
        colorMenu.add(tagOption);
        group.add(tagOption);

        colorMenu.add(getBisulfiteContextMenuItem(group));

        // Base modifications
        mappings.clear();
        mappings.put("base modification", AlignmentTrack.ColorOption.BASE_MODIFICATION);
        mappings.put("base modification (5mC)", AlignmentTrack.ColorOption.BASE_MODIFICATION_5MC);
        mappings.put("base modification (all C)", AlignmentTrack.ColorOption.BASE_MODIFICATION_C);
        colorMenu.addSeparator();
        for (Map.Entry<String, AlignmentTrack.ColorOption> el : mappings.entrySet()) {
            JRadioButtonMenuItem mi = getColorMenuItem(el.getKey(), el.getValue());
            colorMenu.add(mi);
            group.add(mi);
        }

        if (alignmentTrack.getPreferences().getAsBoolean(SMRT_KINETICS_SHOW_OPTIONS)) {
            // Show additional options to help visualize SMRT kinetics data
            mappings.clear();
            mappings.put("SMRT subread IPD", AlignmentTrack.ColorOption.SMRT_SUBREAD_IPD);
            mappings.put("SMRT subread PW", AlignmentTrack.ColorOption.SMRT_SUBREAD_PW);
            mappings.put("SMRT CCS fwd-strand aligned IPD", AlignmentTrack.ColorOption.SMRT_CCS_FWD_IPD);
            mappings.put("SMRT CCS fwd-strand aligned PW", AlignmentTrack.ColorOption.SMRT_CCS_FWD_PW);
            mappings.put("SMRT CCS rev-strand aligned IPD", AlignmentTrack.ColorOption.SMRT_CCS_REV_IPD);
            mappings.put("SMRT CCS rev-strand aligned PW",AlignmentTrack.ColorOption.SMRT_CCS_REV_PW);
            colorMenu.addSeparator();
            for (Map.Entry<String, AlignmentTrack.ColorOption> el : mappings.entrySet()) {
                JRadioButtonMenuItem mi = getColorMenuItem(el.getKey(), el.getValue());
                colorMenu.add(mi);
                group.add(mi);
            }
        }

        add(colorMenu);

    }

    void addShadeAlignmentsMenuItem() {
        JMenu shadeMenu = new JMenu("Shade alignments by");

        for (AlignmentTrack.ShadeAlignmentsOption option : AlignmentTrack.ShadeAlignmentsOption.values()) {
            JRadioButtonMenuItem mi = new JRadioButtonMenuItem(option.label);
            mi.setSelected(renderOptions.getShadeAlignmentsOption() == option);
            mi.addActionListener(aEvt -> {
                alignmentTrack.setShadeAlignmentsOptions(option);
                alignmentTrack.repaint();
            });

            shadeMenu.add(mi);
        }
        add(shadeMenu);
    }

    void addPackMenuItem() {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Re-pack alignments");
        item.addActionListener(aEvt -> UIUtilities.invokeOnEventThread(() -> {
            IGV.getInstance().packAlignmentTracks();
            alignmentTrack.repaint();
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
            final double location = frame.getChromosomePosition(me);

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
            alignmentTrack.setViewAsPairs(viewAsPairs);
        });
        item.setEnabled(dataManager.isPairedEnd());
        add(item);
    }

    void addGoToMate(final TrackClickEvent te, Alignment clickedAlignment) {
        JMenuItem item = new JMenuItem("Go to mate");
        addActionIfMatesAreMapped(te, clickedAlignment, item, this::gotoMate);
        add(item);
    }

    void showMateRegion(final TrackClickEvent te, Alignment clickedAlignment) {
        JMenuItem item = new JMenuItem("View mate region in split screen");
        addActionIfMatesAreMapped(te, clickedAlignment, item, this::splitScreenMate);
        add(item);
    }

    private void addActionIfMatesAreMapped(final TrackClickEvent te, final Alignment clickedAlignment, final JMenuItem item,
                                           final BiConsumer<ReferenceFrame, Alignment> action) {
        final ReferenceFrame frame = te.getFrame();
        if (frame == null) {
            item.setEnabled(false);
        } else {
            final Alignment alignment = clickedAlignment.getSpecificAlignment(te.getChromosomePosition()); // handle PairedAlignments and the like
            item.addActionListener(aEvt -> action.accept(frame, alignment));
            if (alignment == null || !alignment.isPaired() || !alignment.getMate().isMapped()) {
                item.setEnabled(false);
            }
        }
    }

    void addClearSelectionsMenuItem() {
        // Change track height by attribute
        JMenuItem item = new JMenuItem("Clear selections");
        item.addActionListener(aEvt -> {
            alignmentTrack.getSelectedReadNames().clear();
            alignmentTrack.repaint();
        });
        add(item);
    }

    JMenuItem addShowAllBasesMenuItem() {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Show all bases");

        if (renderOptions.getColorOption() == AlignmentTrack.ColorOption.BISULFITE || renderOptions.getColorOption() == AlignmentTrack.ColorOption.NOMESEQ) {
            //    item.setEnabled(false);
        } else {
            item.setSelected(renderOptions.isShowAllBases());
        }
        item.addActionListener(aEvt -> {
            renderOptions.setShowAllBases(item.isSelected());
            alignmentTrack.repaint();
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
            alignmentTrack.repaint();
        });
        add(item);
    }

    JMenuItem addShowMismatchesMenuItem() {
        // Change track height by attribute
        final JMenuItem item = new JCheckBoxMenuItem("Show mismatched bases");
        item.setSelected(renderOptions.isShowMismatches());
        item.addActionListener(aEvt -> {
            renderOptions.setShowMismatches(item.isSelected());
            alignmentTrack.repaint();
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
                alignmentTrack.repaint();
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
            alignmentTrack.repaint();
        }));
        add(item);
    }

    void addShowItems() {

        final CoverageTrack coverageTrack = alignmentTrack.getCoverageTrack();
        if (coverageTrack != null) {
            final JMenuItem item = new JCheckBoxMenuItem("Show Coverage Track");
            item.setSelected(coverageTrack.isVisible());
            item.setEnabled(!coverageTrack.isRemoved());
            item.addActionListener(aEvt -> {
                coverageTrack.setVisible(item.isSelected());
                IGV.getInstance().repaint(Arrays.asList(coverageTrack));

            });
            add(item);
        }

        final SpliceJunctionTrack spliceJunctionTrack = alignmentTrack.getSpliceJunctionTrack();
        if (spliceJunctionTrack != null) {
            final JMenuItem item = new JCheckBoxMenuItem("Show Splice Junction Track");
            item.setSelected(spliceJunctionTrack.isVisible());
            item.setEnabled(!spliceJunctionTrack.isRemoved());
            item.addActionListener(aEvt -> {
                alignmentTrack.setVisible(item.isSelected());
                IGV.getInstance().repaint(Arrays.asList(spliceJunctionTrack));

            });
            add(item);
        }

        final JMenuItem alignmentItem = new JCheckBoxMenuItem("Show Alignment Track");
        alignmentItem.setSelected(true);
        alignmentItem.addActionListener(e -> {
            alignmentTrack.setVisible(alignmentItem.isSelected());
            IGV.getInstance().repaint(Arrays.asList(alignmentTrack));
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
        final Alignment alignment = getSpecificAlignment(te);
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

        final Alignment alignment = getSpecificAlignment(te);
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
        final Alignment alignment = getSpecificAlignment(te);
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

        final Alignment alignment = alignmentTrack.getAlignmentAt(te);
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
        supplementalItem.setSelected(alignmentTrack.isLinkedReads() && "READNAME".equals(renderOptions.getLinkByTag()));
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
        item.setSelected(alignmentTrack.isLinkedReadView() && tag != null && tag.equals(renderOptions.getLinkByTag()));
        item.addActionListener(aEvt -> {
            boolean linkedReads = item.isSelected();
            alignmentTrack.setLinkedReadView(linkedReads, tag);
        });
        return item;
    }

    private JCheckBoxMenuItem linkedReadItem(String tag) {
        final JCheckBoxMenuItem item = new JCheckBoxMenuItem("Link by " + tag);
        item.setSelected(!alignmentTrack.isLinkedReadView() && alignmentTrack.isLinkedReads() && tag.equals(renderOptions.getLinkByTag()));
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

    void addThirdGenItems(Alignment clickedAlignment, final TrackClickEvent tce) {

        //Supplementary/chimeric items, only if the read has an SA tag
        addShowChimericRegions(alignmentTrack, tce, clickedAlignment);


        final JMenuItem qcItem = new JCheckBoxMenuItem("Quick consensus mode");
        qcItem.setSelected(renderOptions.isQuickConsensusMode());
        qcItem.addActionListener(aEvt -> {
            renderOptions.setQuickConsensusMode(qcItem.isSelected());
            alignmentTrack.repaint();
        });

        final JMenuItem thresholdItem = new JMenuItem("Small indel threshold...");
        thresholdItem.addActionListener(evt -> UIUtilities.invokeOnEventThread(() -> {
            String sith = MessageUtils.showInputDialog("Small indel threshold: ", String.valueOf(renderOptions.getSmallIndelThreshold()));
            try {
                renderOptions.setSmallIndelThreshold(Integer.parseInt(sith));
                alignmentTrack.repaint();
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
            alignmentTrack.repaint();
        }));

        final JMenuItem imItem = new JCheckBoxMenuItem("Show insertion markers");
        imItem.setSelected(renderOptions.isShowInsertionMarkers());
        imItem.addActionListener(aEvt -> {
            renderOptions.setShowInsertionMarkers(imItem.isSelected());
            alignmentTrack.repaint();
        });

        add(imItem);
        add(qcItem);
        add(item);
        add(thresholdItem);
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
    private void gotoMate(final ReferenceFrame frame, Alignment alignment) {
        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {

                alignmentTrack.setSelectedAlignment(alignment);

                String chr = mate.getChr();
                int start = mate.start - 1;

                // Don't change scale
                double range = frame.getEnd() - frame.getOrigin();
                int newStart = (int) Math.max(0, (start + (alignment.getEnd() - alignment.getStart()) / 2 - range / 2));
                int newEnd = newStart + (int) range;
                frame.jumpTo(chr, newStart, newEnd);
                frame.recordHistory();
            } else {
                MessageUtils.showMessage("Alignment does not have mate, or it is not mapped.");
            }
        }
    }


    /**
     * Split the screen so the current view and mate region are side by side.
     * Need a better name for this method.
     */
    private void splitScreenMate(ReferenceFrame frame, Alignment alignment) {
        if (alignment != null) {
            ReadMate mate = alignment.getMate();
            if (mate != null && mate.isMapped()) {
                alignmentTrack.setSelectedAlignment(alignment);
                addNewLociToFrames(frame, List.of(mate));
            } else {
                MessageUtils.showMessage("Alignment does not have mate, or it is not mapped.");
            }
        }
    }

    private static void addNewLociToFrames(final ReferenceFrame frame, final List<? extends Locatable> toIncludeInSplit) {
        final List<String> newLoci = toIncludeInSplit.stream()
                .map(locatable -> getLocusStringForAlignment(frame, locatable))
                .collect(Collectors.toList());
        List<String> loci = createLociList(frame, newLoci);
        String listName = String.join("   ", loci); // TODO check the trailing "   " was unnecessary
        //Need to sort the frames by position
        GeneList geneList = new GeneList(listName, loci, false);
        geneList.sort(FrameManager.FRAME_COMPARATOR);
        IGV.getInstance().getSession().setCurrentGeneList(geneList);
        IGV.getInstance().resetFrames();
    }

    private static List<String> createLociList(final ReferenceFrame frame, final List<String> lociToAdd) {
        final List<String> loci = new ArrayList<>(FrameManager.getFrames().size() + lociToAdd.size());
        if (FrameManager.isGeneListMode()) {
            for (ReferenceFrame ref : FrameManager.getFrames()) {
                //If the frame-name is a locus, we use it unaltered
                //Don't want to reprocess, easy to get off-by-one
                String name = ref.getName();
                loci.add(Locus.fromString(name) != null ? name : ref.getFormattedLocusString());
            }
        } else {
            loci.add(frame.getFormattedLocusString());
        }
        loci.addAll(lociToAdd);
        return loci;
    }

    private static String getLocusStringForAlignment(final ReferenceFrame frame, final Locatable alignment) {
        int adjustedMateStart = alignment.getStart() - 1;

        // Generate a locus string for the alignment.  Keep the window width (in base pairs) == to the current range
        Range range = frame.getCurrentRange();
        int length = range.getLength();
        int start = Math.max(0, adjustedMateStart - length / 2);
        int end = start + length;
        return alignment.getContig() + ":" + start + "-" + end;
    }


    /**
     * Get the most "specific" alignment at the specified location.  Specificity refers to the smallest alignment
     * in a group that contains the location (i.e. if a group of linked alignments overlap take the smallest one).
     *
     * @param te
     * @return
     */
    Alignment getSpecificAlignment(TrackClickEvent te) {
        Alignment alignment = alignmentTrack.getAlignmentAt(te);
        if (alignment != null) {
            alignment = alignment.getSpecificAlignment(te.getChromosomePosition());
        }
        return alignment;
    }


    /**
     * Link alignments by arbitrary tag, without the extra settings applied to link-read-view
     *
     * @param linkReads
     * @param tag
     */
    private void setLinkByTag(boolean linkReads, String tag) {
        if (alignmentTrack.isLinkedReadView()) {
            alignmentTrack.undoLinkedReadView();
        }
        if (linkReads) {
            renderOptions.setLinkByTag(tag);
            if (renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.NONE) {
                renderOptions.setGroupByOption(AlignmentTrack.GroupOption.LINKED);
            }
        } else {
            renderOptions.setLinkByTag(null);
            if (renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.LINKED) {
                renderOptions.setGroupByOption(AlignmentTrack.GroupOption.NONE);
            }
        }
        renderOptions.setLinkedReads(linkReads);
        dataManager.packAlignments(renderOptions);
        repaint();
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

}
